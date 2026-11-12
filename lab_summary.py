"""Per-run functional synthesis: signed ORA, multiplicity control, term clustering."""
from pathlib import Path
from functools import lru_cache
import copy,json,re
import numpy as np
from scipy.stats import spearmanr
from scipy.cluster.hierarchy import linkage,fcluster
from scipy.spatial.distance import squareform
from lab_core import bh_adjust,clean_json
from lab_explore import enrichment,selected_features,read_catalogs

@lru_cache(maxsize=1)
def references():
 p=Path(__file__).resolve().parent/'summary_references.json'
 if not p.exists():p=p.parent/'data/summary_references.json'
 return json.loads(p.read_text())

CATEGORIES=[
 ('Gene-expression and biosynthetic regulation',r'regulation of.*biosynth|regulation of.*gene expression|regulation of.*transcription|nucleic acid', 'regulation of gene expression or biosynthesis differs across the samples', 'Measure the highlighted regulators and their downstream RNA and protein products in independent matched samples.'),
 ('Cellular homeostasis',r'homeostas|cellular response to stimulus|response to external stimulus', 'cellular maintenance or stimulus-response state differs across the samples', 'Compare the highlighted proteins with tissue composition and a controlled stimulus-response assay.'),
 ('Energy and mitochondrial function',r'mitochond|respirat|oxidative phosphorylation|ATP synth|tricarbox|glycoly', 'energetic demand differs across the samples', 'Measure oxygen consumption, adenosine triphosphate production and the highlighted metabolic enzymes in new matched samples.'),
 ('RNA processing and gene regulation',r'RNA processing|splic|transcription|chromatin|RNA metabolic', 'RNA regulation contributes to the cross-layer pattern', 'Compare isoform-resolved RNA measurements with targeted protein measurements for the highlighted genes.'),
 ('Protein synthesis and turnover',r'ribosom|translat|proteasom|proteol|protein folding|protein catabol|peptid', 'protein production or turnover differs across the samples', 'Measure new protein synthesis and degradation rates, alongside the highlighted transcripts and proteins.'),
 ('Lipid and membrane metabolism',r'lipid|fatty acid|sterol|phospholipid|sphing|membrane biosynth', 'lipid allocation or membrane composition varies with this component', 'Add targeted lipid measurements and examine the highlighted enzymes in independent samples.'),
 ('Stress and immune response',r'immune|defense|oxidative stress|response to stress|peroxisom|antioxid|inflamm', 'stress-response or immune state varies across the samples', 'Measure reactive oxygen species, antioxidant activity and relevant immune markers; include exposure and pathogen metadata when available.'),
 ('Cell communication and adhesion',r'neuron|synap|axon|cell adhesion|cell.cell|signal transduc', 'cell communication or adhesion contributes to this component', 'Test the highlighted genes with cell-localization and functional assays appropriate to the sampled tissue.'),
 ('Transport and cell organization',r'transport|vesicle|cytoskel|localization|organelle|microtub|actin', 'transport or cellular organization differs across the samples', 'Measure subcellular localization or transport activity for the highlighted proteins.'),
 ('Growth and development',r'development|differentiation|cell cycle|prolifer|growth', 'developmental or proliferative state varies across the samples', 'Measure tissue composition and stage-specific markers in new samples, alongside the highlighted genes.'),
 ('Nutrient and small-molecule metabolism',r'amino acid|carbohydrate|purine|nucleotide|organonitrogen|nutrient', 'nutrient use or small-molecule metabolism varies across the samples', 'Quantify the corresponding metabolites and enzyme activities in independently collected matched samples.')]

def category(name):
 for title,pattern,hyp,test in CATEGORIES:
  if re.search(pattern,name,re.I):return title,hyp,test
 if re.search(r'metabolic|catabolic|biosynth',name,re.I):return 'Metabolism and biosynthesis','metabolic or biosynthetic state varies across the samples','Measure the highlighted enzymes together with targeted metabolites in independent samples.'
 return 'Functional annotation profile','the highlighted molecular process varies across the samples','Measure the highlighted genes and the relevant process in new matched samples.'

def kappa(a,b,universe):
 n=len(universe)
 if not n:return 0.
 a,b=set(a)&universe,set(b)&universe
 po=(len(a&b)+len(universe-(a|b)))/n
 pe=(len(a)*len(b)+(n-len(a))*(n-len(b)))/(n*n)
 return (po-pe)/(1-pe) if pe<1 else 0.

def cluster_terms(rows,universe,cutoff=.5):
 if not rows:return []
 d=np.zeros((len(rows),len(rows)))
 for i,a in enumerate(rows):
  for j in range(i):
   b=rows[j];sim=kappa(a['members'],b['members'],universe) if len(set(a['members'])&set(b['members']))>=3 else 0.
   d[i,j]=d[j,i]=1-max(0.,sim)
 labels=fcluster(linkage(squareform(d,checks=False),method='complete'),t=1-cutoff,criterion='distance') if len(rows)>1 else [1]
 groups=[]
 for label in sorted(set(labels)):
  members=[r for r,l in zip(rows,labels) if l==label]
  members.sort(key=lambda r:(r['report_q'],r['p_value'],r['term_in_background'],r['name']))
  lead=members[0];score=float(np.mean([-np.log10(max(r['p_value'],1e-300)) for r in members]))
  groups.append(dict(title=lead['name'],category=category(lead['name'])[0],score=score,
   best_term_q=min(r['q_value'] for r in members),report_q=min(r['report_q'] for r in members),
   status='BH-supported' if any(r['report_q']<.05 for r in members) else 'Candidate',
   views=sorted(set(r['view'] for r in members)),directions=sorted(set(r['direction'] for r in members)),
   genes=sorted(set().union(*(set(r['features']) for r in members))),terms=members))
 return sorted(groups,key=lambda g:(g['report_q'],-g['score'],-len(g['genes'])))

def model_result(result):
 if result.get('mofa'):return result,'MOFA','Factor'
 r=copy.deepcopy(result);pca=r['joint_pca'];z=np.array(pca['scores']);load=np.array(pca['loadings']);off=0;weights=[];inputs=[];r2=[];tot=[]
 for v in r['views']:
  p=v['model_feature_count'];w=load[off:off+p];off+=p;x=np.array(v['pca_matrix'])/np.sqrt(p);x-=x.mean(0);den=np.sum(x*x)
  weights.append(w);inputs.append(x);r2.append([1-np.sum((x-np.outer(z[:,k],w[:,k]))**2)/den for k in range(z.shape[1])]);tot.append(1-np.sum((x-z@w.T)**2)/den)
 r['mofa']=dict(factors=z.shape[1],scores=z,weights=weights,model_input=inputs,variance_explained=r2,total_r2=tot)
 return clean_json(r),'Joint PCA','PC'

def factor_associations(r,k,metadata):
 by={m['sample_id']:m for m in metadata};ids=r['sample_ids'];z=np.array(r['mofa']['scores'])[:,k];out=[]
 traits={'group':r['groups']}
 for key in ['year','province','age','sex','rna_library_size','rna_TMM_factor','ms_plate']:
  values=[by.get(s,{}).get(key) for s in ids]
  if sum(v is not None and v!='' for v in values)>=4:traits[key]=values
 for key,vals in traits.items():
  ok=np.array([v is not None and v!='' for v in vals]);levels=list(dict.fromkeys(v for v in vals if v is not None and v!=''))
  try:x=np.array([float(v) if v is not None and v!='' else np.nan for v in vals]);label=key
  except (TypeError,ValueError):
   if len(levels)!=2:continue
   x=np.array([levels.index(v) if v in levels else np.nan for v in vals]);label=f'{key}: {levels[1]} = 1'
  ok&=np.isfinite(x)
  if ok.sum()<4 or np.std(x[ok])==0:continue
  rho=float(spearmanr(x[ok],z[ok]).statistic)
  if np.isfinite(rho):out.append(dict(trait=label,rho=rho,n=int(ok.sum())))
 return sorted(out,key=lambda x:-abs(x['rho']))

def summarize(result,options):
 r,model,prefix=model_result(result);context=options.get('context',{});meta=context.get('metadata',[]);labels=context.get('gene_labels',{});top=int(options.get('top_n',20))
 if not 10<=top<=60:raise ValueError('Choose 10–60 top features for the summary.')
 scope=r.get('annotation_scope','human');ann=read_catalogs(scope);queries=[];all_rows=[]
 allowed={'Proteomics':['GO_BP','Bee_KEGG'] if scope=='bee' else ['GO_BP','Reactome'], 'RNA':['GO_BP','Bee_KEGG'] if scope=='bee' else ['GO_BP','Reactome'],'Metabolomics':['Metabolite_pathways'] if context.get('annotation_mode')!='human_genes_only' else []}
 for k in range(r['mofa']['factors']):
  for i,v in enumerate(r['views']):
   if r['mofa']['variance_explained'][i][k]<.01:continue
   for catalog in allowed.get(v['name'],[]):
    for direction in ['positive','negative']:
     q=enrichment(r,dict(factor=k,view=v['name'],catalog=catalog,direction=direction,top_n=top),ann)
     queries.append(q)
     for row in q['rows']:row.update(factor=k+1,view=v['name'],catalog=catalog,direction=direction);all_rows.append(row)
 for row,q in zip(all_rows,bh_adjust([row['p_value'] for row in all_rows])):row['report_q']=float(q)
 term_lookup={(cat,t['id']):t for cat,ts in ann['catalogs'].items() for t in ts};refs=references();factors=[];used={x['id']:x for x in refs if x['kind']=='method'}
 for k in range(r['mofa']['factors']):
  selected=selected_features(r,k,min(top,40));pool=[row for row in all_rows if row['factor']==k+1 and row['overlap']>=2 and (row.get('fold_enrichment') or 0)>1 and row['term_in_background']<=.6*row['background_size']]
  pool.sort(key=lambda x:(x['report_q'],x['p_value'],-x['overlap']))
  candidates=[row for row in pool if row['p_value']<=.1][:80]
  # A ranked lead is still useful when a small gene universe has no low-p terms.
  if not candidates:candidates=pool[:8]
  for row in candidates:row['members']=term_lookup[(row['catalog'],row['id'])]['members']
  gene_universe=set().union(*(set(v['selected_features']) for v in r['views'] if v['name'] in ['Proteomics','RNA']))
  # Compounds are clustered separately so they cannot collide with gene identifiers.
  gene_rows=[x for x in candidates if x['view']!='Metabolomics'];compound_rows=[x for x in candidates if x['view']=='Metabolomics']
  groups=cluster_terms(gene_rows,gene_universe)+cluster_terms(compound_rows,set().union(*(set(v['selected_features']) for v in r['views'] if v['name']=='Metabolomics')))
  groups=sorted(groups,key=lambda g:(g['report_q'],-g['score']))[:5]
  hypotheses=[]
  for g in groups[:2]:
   cat,hyp,test=category(g['title']);genes=', '.join(labels.get(gene,gene) for gene in g['genes'][:5]);topic=' '.join(t['name'] for t in g['terms']).lower()
   related=[ref for ref in refs if ref['kind']=='biological' and ref['species']==scope and any(tag.lower() in topic for tag in ref['tags'])][:2]
   for ref in related:used[ref['id']]=ref
   contrib={v['name']:r['mofa']['variance_explained'][i][k] for i,v in enumerate(r['views'])};pr=contrib.get('Proteomics',0);rr=contrib.get('RNA',0)
   statement=f"One hypothesis is that {hyp}, involving the highlighted {g['title'].lower()} genes/features."
   if 'RNA' in g['views'] and rr>=.05 and rr>3*max(pr,.001):
    statement=f"{prefix} {k+1} reconstructs more of the RNA input ({rr:.1%}) than the protein input ({pr:.1%}). The {g['title'].lower()} signature could be buffered between RNA expression and protein accumulation."
    test='Track the highlighted genes at RNA and protein levels across matched time points, and compare the timing and size of their responses.'
   elif 'Proteomics' in g['views'] and pr>=.05 and pr>3*max(rr,.001):
    statement=f"{prefix} {k+1} reconstructs more of the protein input ({pr:.1%}) than the RNA input ({rr:.1%}). Differences in translation or protein turnover could contribute to the {g['title'].lower()} signature."
    test='Measure the highlighted transcripts alongside new protein synthesis and protein half-lives in independently collected matched samples.'
   elif all(v in g['views'] for v in ['Proteomics','RNA']) and pr>=.05 and rr>=.05:
    statement=f"Both protein ({pr:.1%}) and RNA ({rr:.1%}) contribute to {prefix} {k+1}. Test whether the {g['title'].lower()} genes show coordinated regulation across the two layers."
   hypotheses.append(dict(title=cat,text=statement,basis=f"{g['title']}; contributing genes/features: {genes}.",test=test,references=[x['id'] for x in related],status=g['status']))
  # Gene-specific reading is kept separate from term-derived functional themes.
  top_ids={x['feature'] for rows in selected.values() for x in rows[:10]}
  if scope=='bee' and context.get('tissue')=='Head' and 'GeneID_408786' in top_ids:
   ref=next(x for x in refs if x['id']=='bee_dscam');used[ref['id']]=ref
   hypotheses.append(dict(title='Dscam follow-up',text=f'Dscam ranks among the strongest features for {prefix} {k+1}. A follow-up could test whether its expression or splice forms track the measured sample differences.',basis='Gene-specific lead: GeneID_408786 (Dscam).',test='Measure Dscam isoforms by targeted RNA sequencing and its protein abundance in independent head samples; pair the measurements with an appropriate learning assay.',references=['bee_dscam'],status='Gene-specific hypothesis'))
  if not hypotheses:hypotheses=[dict(title='Follow the strongest features',text=f'Test whether the leading genes/features for {prefix} {k+1} reproduce in a new matched sample set.',basis='The current annotation profile does not resolve a specific functional theme.',test='Repeat the paired measurements, check sample composition and recorded technical factors, then repeat functional annotation.',references=[],status='Open functional assignment')]
  factor_rows=[x for x in all_rows if x['factor']==k+1]
  factors.append(dict(index=k+1,label=f'{prefix} {k+1}',contributions=[dict(view=v['name'],r2=r['mofa']['variance_explained'][i][k]) for i,v in enumerate(r['views'])],
   total_significant_terms=sum(x['report_q']<.05 for x in factor_rows),themes=groups,hypotheses=hypotheses[:3],
   strongest_features={v:[dict(x,label=labels.get(x['feature'],x['feature'])) for x in rows[:6]] for v,rows in selected.items()},
   associations=factor_associations(r,k,meta)[:4],queries=[{x:q[x] for x in ['view','catalog','direction','query_mapped','query_selected','universe_size','model_features','tested_terms']} for q in queries if q['factor']==k+1]))
 papers=[]
 for p in context.get('papers',[]):
  if str(p.get('url','')).startswith('https://'):papers.append(p)
 return clean_json(dict(schema_version='1.0',action='summary',dataset=r['dataset'],model=model,annotation_provenance=ann.get('provenance',{}),component_name=prefix,samples=len(r['sample_ids']),views=[v['name'] for v in r['views']],parameters=r['parameters'],
  top_n=top,reconstruction_basis='centered variation' if model=='Joint PCA' or r['parameters']['mofa_center'] else 'uncentered input energy',cluster_kappa=.5,multiple_testing_family='All tested factor/view/catalog/direction term queries in this summary',tested_term_queries=len(all_rows),
  method='Signed top-weight gene sets; one-sided hypergeometric over-representation. Within-query and report-wide Benjamini-Hochberg correction. Terms with at least two query hits are grouped using complete-linkage clustering at kappa >= 0.5, with at least three shared model genes. Cluster score is mean(-log10(raw term p)); it summarizes related terms and is not an additional p value. This is a DAVID-inspired local grouping, using ordinary hypergeometric p values rather than DAVID EASE scores.',
  factors=factors,references=list(used.values()),study_papers=papers))
