'use strict';
// Per-study choices come from the reproducible, label-free defaults review.
function applyStudyDefaults(data){
 const id=data.teaching?.id||'ccle',p=validatedDefaults[id]||{method:'mofa',factors:3,model_features:100,feature_scaling:'unit_variance',max_missing:.3};
 for(const [id,value] of Object.entries({transform:'auto',normalization:'none',standardize:p.feature_scaling,'mofa-center':'true','scale-views':'false','factor-count':p.factors,'feature-count':p.model_features,'max-missing':p.max_missing,method:p.method,'gene-alignment':p.feature_alignment||'matched_genes'}))$(id).value=String(value);
 $('gene-alignment').disabled=!['RNA','Proteomics'].every(n=>data.views.some(v=>v.name===n));
 networkSettings={mode:'factor',top_per_view:6,threshold:.6};pathwaySettings={view:'Proteomics',catalog:'GO_BP',mode:'ora',top_n:20,direction:'absolute'};
 if($('restore-defaults'))$('restore-defaults').onclick=()=>{applyStudyDefaults(dataset);challenge='none';stale();};
}
const summaryBeforeReview=renderStudySummary;
renderStudySummary=function(){summaryBeforeReview();const t=dataset.teaching||{},p=validatedDefaults[t.id||'ccle'];if(!p)return;const bridge=$('public-study-overview').querySelector('.mofa-bridge strong');if(bridge)bridge.textContent=p.method==='pca'?'Joint PCA':'MOFA';
 $('public-study-overview').insertAdjacentHTML('beforeend',`<div class="default-explainer"><b>Recommended settings</b><p>${p.method==='pca'?'Start with joint PCA for these six donors. Select MOFA to compare a factor model.':`MOFA · ${p.factors} factors · ${p.feature_scaling==='pareto'?'Pareto':'unit-variance'} scaling · 100 features per view. Settings were selected from three repeated masked-entry comparisons.`} Protein and RNA use verified common genes in the same order. <a href="sources.html?v=20261112s#gene-alignment">How matching and defaults work ↗</a></p></div>`);
};
const metadataBeforeReview=metadata;
metadata=function(){return metadataBeforeReview()+`<span>${result.gene_alignment?.mode==='matched_genes'?result.gene_alignment.model_pairs.length+' protein–RNA genes aligned by ID and column order':'Independent view features'}</span>`;};
const integrationBeforeReview=renderIntegration;
renderIntegration=function(box){integrationBeforeReview(box);
 const a=result.gene_alignment;if(!a?.model_pairs?.length)return;
 box.querySelector('.integration-intro p').textContent=`${a.model_pairs.length} verified genes are represented in both protein and RNA, in exactly the same order. Sample rows are matched by ID. ${result.views.length>2?'Metabolites retain their chemical feature identities.':''}`;
 const pairs=a.model_pairs;
 box.insertAdjacentHTML('beforeend',`<section class="gene-alignment-panel"><span class="kicker">GENE MATCHING</span><h2>Protein–RNA gene alignment</h2><p>${a.eligible_pairs} shared genes pass quality checks in both layers; ${pairs.length} enter both model matrices. Gene selection combines within-view variance ranks without using group labels.</p><div class="table-wrap"><table><thead><tr><th>Column</th><th>Gene</th><th>Protein feature ID</th><th>RNA feature ID</th></tr></thead><tbody>${pairs.slice(0,6).map((p,i)=>`<tr><td>${i+1}</td><td>${esc(p.gene_symbol)}</td><td>${esc(p.protein_feature_id)}</td><td>${esc(p.rna_feature_id)}</td></tr>`).join('')}</tbody></table></div><p class="note">The complete gene crosswalk is included in the result bundle. Each assay retains its own measured values.</p>${result.mofa?'<label class="factor-control">Inspect matched genes<select id="paired-factor">'+factorOptions()+'</select></label>'+chart('paired-gene-weights','The same genes, side by side','Factor 1 · top 12 matched genes by joint absolute weight. Each view’s weights are divided by that view’s maximum absolute weight. Actual weights are shown on hover.','tall'):''}</section>`);
 if(result.mofa){const m=result.mofa,pi=result.views.findIndex(v=>v.name==='Proteomics'),ri=result.views.findIndex(v=>v.name==='RNA'),pm=new Map(result.views[pi].selected_features.map((f,i)=>[f,i])),rm=new Map(result.views[ri].selected_features.map((f,i)=>[f,i]));
 function drawPaired(){const k=Number($('paired-factor').value);const max=[pi,ri].map(i=>Math.max(1e-12,...m.weights[i].map(w=>Math.abs(w[k]))));
 const rows=pairs.map(p=>({p,w:[m.weights[pi][pm.get(p.protein_feature_id)][k],m.weights[ri][rm.get(p.rna_feature_id)][k]]})).map(r=>({...r,z:r.w.map((w,i)=>w/max[i])})).sort((a,b)=>b.z.reduce((s,x)=>s+Math.abs(x),0)-a.z.reduce((s,x)=>s+Math.abs(x),0)).slice(0,12);
 plot('paired-gene-weights',[{type:'heatmap',x:['Protein','RNA'],y:rows.map(r=>r.p.gene_symbol),z:rows.map(r=>r.z),customdata:rows.map(r=>r.w),text:rows.map(r=>r.z.map(w=>fmt(w,2))),texttemplate:'%{text}',zmin:-1,zmax:1,colorscale:[[0,'#287d9e'],[.5,'#ffffff'],[1,'#bd683c']],colorbar:{title:'Relative<br>weight'},hovertemplate:'%{y}<br>%{x}<br>Actual weight %{customdata:.4f}<br>Relative weight %{z:.3f}<extra></extra>'}],{margin:{l:180,r:90,t:15,b:65},xaxis:{title:'Molecular view'},yaxis:{title:'Matched gene',autorange:'reversed'}});$('paired-gene-weights').parentElement.querySelector('p').textContent='Factor '+(k+1)+' · 12 matched genes with the largest joint relative weights. The row always represents the same gene in both views.';bindPlotSaves();}
 $('paired-factor').value=String(factorIndex);$('paired-factor').onchange=drawPaired;drawPaired();
 }
 bindPlotSaves();
};

// A completed model opens each explorer with its reviewed starting query.
const networkBeforeReview=renderNetworkPanel;
renderNetworkPanel=function(box){networkBeforeReview(box);if(result.mofa&&!networkResult&&!busy)queueMicrotask(()=>{if(current==='network'&&!busy)$('update-network')?.click();});};
const pathwaysBeforeReview=renderPathwaysPanel;
renderPathwaysPanel=function(box){pathwaysBeforeReview(box);if(result.mofa&&!enrichmentResult&&!busy)queueMicrotask(()=>{if(current==='pathways'&&!busy)$('update-pathways')?.click();});};
const endRequestBeforeReview=endRequest;
endRequest=function(){endRequestBeforeReview();queueMicrotask(()=>{if(busy||!result?.mofa)return;if(current==='network'&&!networkResult)$('update-network')?.click();if(current==='pathways'&&!enrichmentResult)$('update-pathways')?.click();});};
// Open each analysis chapter at its heading; retain position for tissue changes.
const chapterBeforeScroll=setChapter;
setChapter=function(id){const previous=current;chapterBeforeScroll(id);if(current!==previous)window.scrollTo({top:0,left:0,behavior:'instant'});};
