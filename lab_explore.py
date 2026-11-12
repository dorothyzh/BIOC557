"""Source-backed exploratory interpretation and real K sensitivity for BIOC557."""
from pathlib import Path
from dataclasses import asdict, replace
import json
from functools import lru_cache
import numpy as np
from scipy.stats import hypergeom, spearmanr
from scipy.optimize import linear_sum_assignment
from lab_core import analyze, _fit_mofa, Parameters, bh_adjust, clean_json

@lru_cache(maxsize=2)
def read_catalogs(scope='human'):
    if scope not in ['human','bee']:raise ValueError('Unknown annotation species.')
    name='bee_enrichment_annotations.json' if scope=='bee' else 'enrichment_annotations.json'
    path=Path(__file__).resolve().parent/name
    if not path.exists():path=path.parent/'data'/name
    return json.loads(path.read_text())

def compare_counts(data, parameters, counts=(2,3,4,5,6)):
    counts=sorted(set(int(k) for k in counts))
    if not counts or len(counts)>5 or min(counts)<2 or max(counts)>6:
        raise ValueError('Compare2–6 factors, at most5 models.')
    p=dict(parameters,method='pca')
    p['factors']=max(counts)
    base=analyze(data,p)
    params=Parameters(**base['parameters'])
    matrices=[np.array(v['processed_matrix'],dtype=float) for v in base['views']]
    features=[v['selected_features'] for v in base['views']]
    names=[v['name'] for v in base['views']]
    models=[]
    for k in counts:
        fitted=_fit_mofa(matrices,features,base['sample_ids'],names,replace(params,factors=k))
        models.append({'requested_k':k,'factors':fitted['factors'],'total_r2':fitted['total_r2'],
                       'factor_r2':fitted['variance_explained'],'scores':fitted['scores']})
    ref_k=int(parameters.get('factors',3))
    reference=next((m for m in models if m['requested_k']==ref_k),models[0])
    for model in models:
        correlations=np.zeros((reference['factors'],model['factors']))
        for i in range(reference['factors']):
            for j in range(model['factors']):
                a,b=np.array(reference['scores'])[:,i],np.array(model['scores'])[:,j]
                r=spearmanr(a,b).statistic if np.std(a)>1e-12 and np.std(b)>1e-12 else 0
                correlations[i,j]=abs(r) if np.isfinite(r) else 0
        rows,cols=linear_sum_assignment(-correlations)
        model['matched_to_reference']=[{'reference_factor':int(i+1),'factor':int(j+1),'absolute_score_spearman':correlations[i,j]} for i,j in zip(rows,cols)]
    return clean_json({'models':models,'view_names':names,'reference_k':reference['requested_k'],
        'parameters':dict(parameters),'baseline':'centered variation' if params.mofa_center else 'uncentered input energy',
        'note':'Actual refits with identical features, samples, preprocessing and seed. In-sample reconstruction is not predictive performance or a rule for optimal K. Factor identities are matched by absolute score correlation, not factor number.'})

def selected_features(result,factor=0,top_n=20,direction='absolute'):
    if 'mofa' not in result:raise ValueError('Fit MOFA first.')
    if not 0<=factor<result['mofa']['factors']:raise ValueError('Select an existing factor.')
    if not 3<=top_n<=80:raise ValueError('Select3–80 top features.')
    out={}
    for i,v in enumerate(result['views']):
        w=np.array(result['mofa']['weights'][i],dtype=float)[:,factor]
        if direction=='positive':eligible=np.flatnonzero(w>0);order=eligible[np.argsort(-w[eligible],kind='stable')]
        elif direction=='negative':eligible=np.flatnonzero(w<0);order=eligible[np.argsort(w[eligible],kind='stable')]
        elif direction=='absolute':order=np.argsort(-np.abs(w),kind='stable')
        else:raise ValueError('Choose absolute, positive or negative weights.')
        out[v['name']]=[{'feature':v['selected_features'][j],'weight':float(w[j]),'column':int(j),'view':v['name']} for j in order[:top_n]]
    return out

def enrichment(result, options, annotations):
    factor=int(options.get('factor',0));top_n=int(options.get('top_n',20));direction=options.get('direction','absolute')
    view_name=options.get('view','RNA');catalog=options.get('catalog','GO_BP');mode=options.get('mode','ora')
    view=next((v for v in result['views'] if v['name']==view_name),None)
    if view is None:raise ValueError('Unknown molecular view.')
    if catalog not in annotations['catalogs']:raise ValueError('Unknown annotation catalog.')
    if (view_name=='Metabolomics') != (catalog=='Metabolite_pathways'):
        raise ValueError('GO/Reactome gene sets use genes; metabolites use mapped compound pathways.')
    terms=annotations['catalogs'][catalog]
    annotated=set().union(*(set(t['members']) for t in terms)) if terms else set()
    # ORA background is the selected model features actually eligible for ranking.
    universe=set(view['selected_features']) & annotated
    picks=selected_features(result,factor,top_n,direction)[view_name]
    query={x['feature'] for x in picks}&universe
    min_size=2 if view_name=='Metabolomics' else 3
    rows=[]
    if mode=='coverage':
        if view_name!='Metabolomics':raise ValueError('Coverage view is for mapped metabolites.')
        available=set(view.get('retained_feature_ids',view['selected_features']))&annotated
        for t in terms:
            hit=sorted(set(t['members'])&available)
            if hit:rows.append({'id':t['id'],'name':t['name'],'overlap':len(hit),'features':hit,'source_url':t.get('source_url'),'coverage_denominator':len(available)})
        rows.sort(key=lambda x:(-x['overlap'],x['name']))
        return {'mode':'coverage','catalog':catalog,'view':view_name,'rows':rows,'measured_mapped':len(available),'measured_total':view['retained_features'],
                'note':'Coverage counts mapped measured compounds. It is not enrichment, pathway activity or a hypothesis test.'}
    for term in terms:
        members=set(term['members'])&universe
        if len(members)<min_size:continue
        hit=sorted(query&members);N=len(universe);M=len(members);n=len(query);k=len(hit)
        p=float(hypergeom.sf(k-1,N,M,n)) if n else 1.
        rows.append({'id':term['id'],'name':term['name'],'overlap':k,'query_size':n,'background_size':N,'term_in_background':M,
                     'p_value':p,'features':hit,'fold_enrichment':(k/n)/(M/N) if n and N and M else None,'gene_ratio':k/n if n else 0.,'source_url':term.get('source_url')})
    q=bh_adjust([r['p_value'] for r in rows])
    for row,value in zip(rows,q):row['q_value']=float(value)
    rows.sort(key=lambda x:(x['q_value'],x['p_value'],-x['overlap']))
    return clean_json({'mode':'ora','view':view_name,'catalog':catalog,'factor':factor+1,'direction':direction,'top_n':top_n,
        'universe_size':len(universe),'model_features':len(view['selected_features']),'query_mapped':len(query),'query_selected':len(picks),
        'tested_terms':len(rows),'significant_terms':sum(r['q_value']<.05 for r in rows),'rows':rows,
        'note':'One-sided hypergeometric over-representation; BH across every tested catalog term for this factor/view/query. Background is model features with pathway annotations. Exploratory annotation of selected weights; not a pathway-activity test.'})

def network(result,options):
    factor=int(options.get('factor',0));top=int(options.get('top_per_view',6));mode=options.get('mode','factor');threshold=float(options.get('threshold',.6))
    if not 3<=top<=10:raise ValueError('Choose3–10 features per view.')
    if mode not in ['factor','correlation']:raise ValueError('Choose factor or correlation network.')
    if not .3<=threshold<=.95:raise ValueError('Correlation threshold must be0.3–0.95.')
    selected=selected_features(result,factor,top)
    nodes=[];values=[]
    for i,v in enumerate(result['views']):
        x=np.array(result['mofa']['model_input'][i],dtype=float)
        for row in selected[v['name']]:
            nodes.append(dict(row,id=v['name']+'::'+row['feature'],factor_r2=result['mofa']['variance_explained'][i][factor]))
            values.append(x[:,row['column']])
    edges=[]
    if mode=='factor':
        for node in nodes:edges.append({'source':'Factor '+str(factor+1),'target':node['id'],'value':node['weight'],'kind':'model_weight'})
    else:
        for i,a in enumerate(nodes):
            for j in range(i+1,len(nodes)):
                b=nodes[j]
                if a['view']==b['view']:continue
                x,y=values[i],values[j];keep=np.isfinite(x)&np.isfinite(y);n=int(keep.sum())
                if n<10 or np.std(x[keep])<1e-12 or np.std(y[keep])<1e-12:continue
                rho=float(spearmanr(x[keep],y[keep]).statistic)
                if abs(rho)>=threshold:edges.append({'source':a['id'],'target':b['id'],'value':rho,'n':n,'kind':'observed_spearman'})
        edges=sorted(edges,key=lambda e:-abs(e['value']))[:100]
    return clean_json({'factor':factor+1,'mode':mode,'threshold':threshold,'nodes':nodes,'edges':edges,
        'note':'Factor links are learned model weights; correlation links are observed pairwise-complete Spearman correlations among selected cross-view features (at least10 samples), not physical interactions or causal links.'})

def explore(result,options,annotations=None):
    action=options.get('action','network')
    if action=='summary':
        from lab_summary import summarize
        return summarize(result,options)
    if action=='network':return network(result,options)
    if action=='enrichment':return enrichment(result,options,annotations if annotations is not None else read_catalogs(result.get('annotation_scope','human')))
    raise ValueError('Unknown exploration action.')

def compare_json(dataset_json,parameters_json,counts_json='[2,3,4,5,6]'):
    return json.dumps(compare_counts(json.loads(dataset_json),json.loads(parameters_json),json.loads(counts_json)),allow_nan=False)

def explore_json(result_json,options_json):
    return json.dumps(explore(json.loads(result_json),json.loads(options_json)),allow_nan=False)
