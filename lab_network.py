"""Integrated feature/annotation graph built from the current fitted measurements."""
import numpy as np
from scipy.stats import spearmanr, hypergeom
from lab_core import bh_adjust, clean_json
from lab_explore import read_catalogs


def graph_features(result, factor, top):
    nodes=[]; values=[]
    for vi,view in enumerate(result['views']):
        w=np.asarray(result['mofa']['weights'][vi],float)
        r2=np.maximum(0,np.asarray(result['mofa']['variance_explained'][vi],float))
        if factor>=0:
            importance=np.abs(w[:,factor]); strongest=np.full(w.shape[0],factor,dtype=int)
        else:
            norm=np.maximum(1e-12,np.max(np.abs(w),axis=0))
            weighted=np.abs(w)/norm*np.sqrt(r2)[None,:]
            importance=np.sum(weighted,axis=1);strongest=np.argmax(weighted,axis=1)
        scale=max(1e-12,float(np.max(importance)))
        x=np.asarray(view['processed_matrix'],float)
        for j in np.argsort(-importance,kind='stable')[:top]:
            k=int(strongest[j]);feature=view['selected_features'][j]
            nodes.append(dict(id=view['name']+'::'+feature,kind='feature',feature=feature,label=feature,
                view=view['name'],column=int(j),weight=float(w[j,k]),factor=k+1,
                factor_r2=float(r2[k]),importance=float(importance[j]/scale)))
            values.append(x[:,j])
    return nodes,values


def functional_nodes(result, features, limit, annotations=None):
    if not limit or result.get('annotation_scope')=='none':return [],[],[],{}
    ann=annotations if annotations is not None else read_catalogs(result.get('annotation_scope','human'))
    catalogs=ann['catalogs'];tests=[];mapped={};by_view={v['name']:v for v in result['views']}
    for name,view in by_view.items():
        chosen={n['feature'] for n in features if n['view']==name}
        names=['Metabolite_pathways'] if name=='Metabolomics' else [x for x in ['GO_BP','Reactome','Bee_KEGG'] if x in catalogs]
        for catalog in names:
            terms=catalogs.get(catalog,[])
            annotated=set().union(*(set(t['members']) for t in terms)) if terms else set()
            universe=set(view['selected_features'])&annotated;query=chosen&universe
            mapped[name]=sorted(set(mapped.get(name,[]))|query)
            N=len(universe);n=len(query)
            for term in terms:
                members=set(term['members'])&universe;M=len(members)
                if M<(2 if name=='Metabolomics' else 3):continue
                hits=sorted(query&members);k=len(hits)
                tests.append(dict(id=term['id'],name=term['name'],catalog=catalog,view=name,
                    p_value=float(hypergeom.sf(k-1,N,M,n)) if n else 1.,overlap=k,
                    background_size=N,query_size=n,term_in_background=M,features=hits,
                    source_url=term.get('source_url','')))
    for row,q in zip(tests,bh_adjust([row['p_value'] for row in tests])):row['q_value']=float(q)
    groups={}
    for row in tests:
        if row['overlap']<2 or (row['catalog']!='Metabolite_pathways' and row['term_in_background']>.6*row['background_size']):continue
        key=row['catalog']+'::'+row['id'];groups.setdefault(key,[]).append(row)
    candidates=[]
    for key,rows in groups.items():
        term=next(t for t in catalogs[rows[0]['catalog']] if t['id']==rows[0]['id'])
        memberships=sorted(n['id'] for n in features if n['feature'] in set(term['members']) and (n['view']=='Metabolomics')==(rows[0]['catalog']=='Metabolite_pathways'))
        best=min(rows,key=lambda r:(r['q_value'],r['p_value'],r['term_in_background'],r['name']))
        candidates.append(dict(id='term::'+key,kind='function',feature=best['id'],label=best['name'],
            catalog=best['catalog'],view='Function',q_value=min(r['q_value'] for r in rows),
            p_value=min(r['p_value'] for r in rows),status='BH-supported' if any(r['q_value']<.05 for r in rows) else 'Annotation match',
            memberships=memberships,views=sorted({n['view'] for n in features if n['id'] in memberships}),
            source_url=best['source_url'],queries=rows,overlap=len(memberships)))
    candidates.sort(key=lambda n:(n['q_value'],n['p_value'],-len(n['views']),-n['overlap'],n['label']))
    picked=[]
    def add(node):
        if len(picked)>=limit or node in picked:return
        for old in picked:
            a,b=set(node['memberships']),set(old['memberships'])
            if node['catalog']==old['catalog'] and len(a&b)/max(1,len(a|b))>=.8:return
        picked.append(node)
    # Represent each available annotation catalog, including compound pathways.
    for catalog in dict.fromkeys(n['catalog'] for n in candidates):
        add(next(n for n in candidates if n['catalog']==catalog))
    for node in candidates:add(node)
    edges=[dict(source=feature,target=node['id'],kind='annotation',value=1.,catalog=node['catalog'])
        for node in picked for feature in node['memberships']]
    return picked,edges,tests,mapped


def graph_layout(nodes,edges,view_names):
    n=len(nodes);anchors=np.zeros((n,2));positions=np.zeros((n,2));kinds=np.array([node['kind']=='function' for node in nodes])
    centers={}
    for i,name in enumerate(view_names):
        angle=np.pi*.8-i*2*np.pi/max(2,len(view_names))
        centers[name]=np.array([1.45*np.cos(angle),1.25*np.sin(angle)-.35])
        indices=[j for j,node in enumerate(nodes) if node.get('view')==name]
        for k,j in enumerate(indices):
            theta=(k+1)*2.399963229728653+i*.4;radius=.22+.42*np.sqrt((k+.5)/len(indices))
            anchors[j]=centers[name];positions[j]=anchors[j]+radius*np.array([np.cos(theta),np.sin(theta)])
    terms=np.flatnonzero(kinds)
    for k,j in enumerate(terms):
        columns=min(4,len(terms));row=k//columns;col=k%columns
        anchors[j]=positions[j]=[(col-(columns-1)/2)*1.02,1.65+row*.66]
    lookup={node['id']:i for i,node in enumerate(nodes)}
    for iteration in range(180):
        delta=positions[:,None,:]-positions[None,:,:];dist2=np.sum(delta*delta,axis=2)+.04
        repulsion=.004*np.sum(delta/(dist2[:,:,None]**1.5),axis=1)
        force=repulsion+(anchors-positions)*np.where(kinds[:,None],.24,.055)
        for edge in edges:
            i,j=lookup[edge['source']],lookup[edge['target']]
            d=positions[j]-positions[i];length=max(1e-8,float(np.linalg.norm(d)))
            strength=.005 if edge['kind']=='annotation' else .008
            spring=strength*(length-.55)*d/length;force[i]+=spring;force[j]-=spring
        positions+=np.clip(force,-.12,.12)*.75
    for node,pos in zip(nodes,positions):node['x']=float(pos[0]);node['y']=float(pos[1])


def integrated_network(result, options, annotations=None):
    model='MOFA';component='Factor'
    if not result.get('mofa'):
        from lab_summary import model_result
        result,model,component=model_result(result)
    factor=int(options.get('factor',-1));top=int(options.get('top_per_view',20));threshold=float(options.get('threshold',.6))
    q_max=float(options.get('q_max',.05));limit=int(options.get('term_limit',8));within=bool(options.get('include_within',False));identity=bool(options.get('include_identity',True))
    if not -1<=factor<result['mofa']['factors']:raise ValueError('Choose all factors or an existing factor.')
    if not 4<=top<=40 or not 0<=limit<=12:raise ValueError('Use 4–40 features per view and 0–12 function nodes.')
    if not .3<=threshold<=.95 or not 0<q_max<=1:raise ValueError('Choose a correlation threshold and BH q cutoff.')
    features,values=graph_features(result,factor,top);tests=[]
    for i,a in enumerate(features):
        for j in range(i+1,len(features)):
            b=features[j]
            if not within and a['view']==b['view']:continue
            x,y=values[i],values[j];ok=np.isfinite(x)&np.isfinite(y);n=int(ok.sum())
            if n<10 or np.std(x[ok])<1e-12 or np.std(y[ok])<1e-12:continue
            test=spearmanr(x[ok],y[ok]);tests.append(dict(source=a['id'],target=b['id'],value=float(test.statistic),p_value=float(test.pvalue),n=n,
                kind='correlation',view_pair=' / '.join(sorted([a['view'],b['view']]))))
    for edge,q in zip(tests,bh_adjust([e['p_value'] for e in tests])):edge['q_value']=float(q)
    eligible=sorted([e for e in tests if abs(e['value'])>=threshold and e['q_value']<=q_max],key=lambda e:(-abs(e['value']),e['q_value'],e['source'],e['target']))
    pairs=list(dict.fromkeys(e['view_pair'] for e in eligible));budget=160
    chosen=[]
    for pair in pairs:chosen.extend([e for e in eligible if e['view_pair']==pair][:max(1,budget//max(1,len(pairs)))])
    for edge in eligible:
        if len(chosen)>=budget:break
        if edge not in chosen:chosen.append(edge)
    edges=chosen[:budget];ids={n['id'] for n in features}
    if identity:
        for pair in result.get('gene_alignment',{}).get('model_pairs',[]):
            source='Proteomics::'+pair['protein_feature_id'];target='RNA::'+pair['rna_feature_id']
            if source in ids and target in ids:edges.append(dict(source=source,target=target,kind='matched_gene',value=1.,gene=pair.get('gene_symbol',pair.get('gene_id',''))))
    functions,membership,term_tests,mapped=functional_nodes(result,features,limit,annotations)
    edges+=membership;nodes=features+functions;graph_layout(nodes,edges,[v['name'] for v in result['views']])
    return clean_json(dict(mode='integrated',dataset=result['dataset'],samples=len(result['sample_ids']),parameters=result['parameters'],factor=factor+1,focus=('All factors' if component=='Factor' else 'All components') if factor<0 else component+' '+str(factor+1),model=model,component_name=component,nodes=nodes,edges=edges,
        options=dict(factor=factor,top_per_view=top,threshold=threshold,q_max=q_max,term_limit=limit,include_within=within,include_identity=identity),
        stats=dict(features=len(features),functions=len(functions),correlation_tests=len(tests),eligible_correlations=len(eligible),displayed_correlations=min(budget,len(chosen)),
            matched_gene_links=sum(e['kind']=='matched_gene' for e in edges),annotation_links=len(membership),functional_tests=len(term_tests),mapped_features={k:len(v) for k,v in mapped.items()}),
        correlation_tests=tests,functional_tests=term_tests,
        methods=dict(correlations='Pairwise observed Spearman correlation, at least 10 pairs. BH across all eligible feature-pair tests in this graph before the rho/q filters. At most 160 links, balanced across view pairs.',
            functions='Hypergeometric over-representation of displayed features against annotated model features in each view. BH across all tested view/catalog terms in this graph. Function nodes show their lowest query q; membership edges show catalog annotations.',
            selection='All-factor ranking sums absolute weights normalized within factor/view and weighted by sqrt(view reconstruction R2). A single-factor view ranks absolute weights.',
            identity='Matched protein/RNA gene identifiers from the verified model crosswalk.')))
