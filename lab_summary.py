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

def reading_for_factor(groups,selected,labels,curated,hypotheses):
 """Curated contextual reading, separate from the single retained interpretation."""
 categories={g['category'] for g in groups};features={x for v in selected.values() for row in v for x in (row['feature'],labels.get(row['feature'],row['feature']))}
 support={p['id']:p['support'] for h in hypotheses for p in h['paper_support']};rows=[]
 for ref in curated:
  spec=ref.get('reading',{});matched=sorted(categories&set(spec.get('categories',[])));genes=sorted(features&set(spec.get('features',[])))
  if not (matched or genes or spec.get('study_wide') or ref['id'] in support):continue
  rows.append(dict(id=ref['id'],relevance=support.get(ref['id'],spec.get('relevance','')),
   role='Interpretation support' if ref['id'] in support else 'Functional context',matched_themes=matched,matched_features=genes,
   priority=100 if ref['id'] in support else spec.get('priority',0)+10*bool(genes)+len(matched)))
 rows.sort(key=lambda x:(-x['priority'],x['id']))
 return [{k:v for k,v in row.items() if k!='priority'} for row in rows[:8]]

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
 term_lookup={(cat,t['id']):t for cat,ts in ann['catalogs'].items() for t in ts};refs=references();factors=[]
 study_key='bee' if scope=='bee' else 'islets' if 'islet' in str(context.get('id','')).lower() or context.get('annotation_mode')=='human_genes_only' else 'ccle'
 curated=[ref for ref in refs if study_key in ref.get('reading',{}).get('studies',[]) and (not ref['reading'].get('tissues') or context.get('tissue') in ref['reading']['tissues'])]
 used={x['id']:x for x in refs if x['kind']=='method' and (not x.get('applies_to') or study_key in x['applies_to'])}
 for ref in curated:used[ref['id']]=ref
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
  # One reviewed interpretation at most. No keyword-to-paper or single-gene fallback.
  hypotheses=[]
  contributions={v['name']:r['mofa']['variance_explained'][i][k] for i,v in enumerate(r['views'])}
  # Reviewed support is tied to specific RNA annotations and tissue, not category keywords.
  reviewed_rna_terms={'GO:0009889','GO:0031326','GO:0010468','GO:0010556','GO:0010467','GO:0009059','GO:0044249','GO:0016070','GO:0051252','GO:0006351','GO:0032774','GO:0006355','GO:2001141','GO:0006357','GO:0006366'}
  support_reasons={
   'bee_crop':'The source study identifies biosynthesis and protein-synthesis modules in matched honeybee RNA and protein measurements.',
   'bee_brain_grn':'Independent honeybee brain data identify coordinated gene-regulatory programs across cell types and behavioral maturation.',
   'bee_brain_cells':'Single-cell brain profiles resolve cell-type and caste-associated RNA expression, providing context for cellular composition in bulk head samples.'}
  centered=model=='Joint PCA' or r['parameters']['mofa_center']
  if centered and len(r['sample_ids'])>=16 and scope=='bee' and context.get('tissue')=='Head':
   for g in groups:
    eligible=[t for t in g['terms'] if t['id'] in reviewed_rna_terms and t['view']=='RNA' and t['report_q']<.05 and t['overlap']>=3 and contributions.get('RNA',0)>=.10]
    ids=('bee_crop','bee_brain_grn','bee_brain_cells')
    if not eligible:continue
    eligible.sort(key=lambda t:(t['report_q'],t['p_value'],-t['overlap']))
    lead=eligible[0];view=lead['view'];lead_genes=lead['features'];cat=g['category']
    statement=f"{prefix} {k+1} captures an RNA expression pattern enriched for {lead['name'].lower()}."
    test='Check whether the same regulatory-gene pattern recurs in independent head pools, with crop and collection year balanced.'
    associations=factor_associations(r,k,meta)
    year=next((x for x in associations if x['trait']=='year' and abs(x['rho'])>=.5),None)
    if year:statement+=f" Its scores are associated with collection year (Spearman rho {year['rho']:.2f})."
    paper_support=[]
    for id in ids:
     ref=next((x for x in refs if x['id']==id),None)
     if ref:used[id]=ref;paper_support.append(dict(id=id,support=support_reasons[id]))
    if len(paper_support)<2:continue
    genes=', '.join(labels.get(gene,gene) for gene in lead_genes[:8])
    hypotheses=[dict(title=cat,text=statement,basis=f"{lead['name']}; {view} contribution {contributions[view]:.1%}; {lead['overlap']} selected genes; report-wide BH q {lead['report_q']:.3g}. Genes: {genes}.",
     test=test,references=list(ids),paper_support=paper_support,status='Enrichment- and literature-supported',evidence_term=lead['id'],evidence_view=view,report_q=lead['report_q'])]
    break
  factor_rows=[x for x in all_rows if x['factor']==k+1]
  support_ids={id for h in hypotheses for id in h['references']}
  factor_reading=reading_for_factor(groups,selected,labels,curated+[ref for ref in refs if ref['id'] in support_ids and ref not in curated],hypotheses)
  factors.append(dict(index=k+1,label=f'{prefix} {k+1}',contributions=[dict(view=v['name'],r2=r['mofa']['variance_explained'][i][k]) for i,v in enumerate(r['views'])],
   total_significant_terms=sum(x['report_q']<.05 for x in factor_rows),themes=groups,hypotheses=hypotheses[:1],related_reading=factor_reading,
   strongest_features={v:[dict(x,label=labels.get(x['feature'],x['feature'])) for x in rows[:6]] for v,rows in selected.items()},
   associations=factor_associations(r,k,meta)[:4],queries=[{x:q[x] for x in ['view','catalog','direction','query_mapped','query_selected','universe_size','model_features','tested_terms']} for q in queries if q['factor']==k+1]))
 papers=[]
 for p in context.get('papers',[]):
  if str(p.get('url','')).startswith('https://'):
   doi=str(p.get('doi','')).lower() or str(p['url']).lower().replace('https://doi.org/','')
   canonical=next((ref for ref in refs if ref['doi'].lower()==doi or ref['title'].lower()==p.get('title','').lower()),None)
   papers.append(dict(p,**{key:canonical[key] for key in ['title','journal','year','year_display','url','doi'] if canonical and key in canonical}))
 return clean_json(dict(schema_version='1.0',action='summary',dataset=r['dataset'],model=model,annotation_provenance=ann.get('provenance',{}),component_name=prefix,samples=len(r['sample_ids']),views=[v['name'] for v in r['views']],parameters=r['parameters'],
  interpretation_policy='At most one interpretation; report-wide BH q < 0.05, at least three selected supporting genes, at least 10% reconstruction in the supporting view, at least 16 samples, centered inputs and two explicitly reviewed relevant primary studies. The current reviewed interpretation covers specific RNA-regulatory annotations in honeybee head samples; other contexts show functional results. No single-gene or keyword-only citation fallback.',top_n=top,reconstruction_basis='centered variation' if model=='Joint PCA' or r['parameters']['mofa_center'] else 'uncentered input energy',cluster_kappa=.5,multiple_testing_family='All tested factor/view/catalog/direction term queries in this summary',tested_term_queries=len(all_rows),
  method='Signed top-weight gene sets; one-sided hypergeometric over-representation. Within-query and report-wide Benjamini-Hochberg correction. Terms with at least two query hits are grouped using complete-linkage clustering at kappa >= 0.5, with at least three shared model genes. Cluster score is mean(-log10(raw term p)); it summarizes related terms and is not an additional p value. This is a DAVID-inspired local grouping, using ordinary hypergeometric p values rather than DAVID EASE scores.',
  related_reading_policy='Primary research is curated by study, tissue, functional theme and measured features. Interpretation support is labelled separately from functional context.',
  reading_library=[dict(id=ref['id'],relevance=ref['reading']['relevance']) for ref in curated],
  factors=factors,references=list(used.values()),study_papers=papers))
