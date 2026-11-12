'use strict';
// One tissue per analysis, with a species-specific annotation universe.
const beeEntry={id:'bee-own-head-2025',path:'data/demos/bee-own-head-2025.json',tag:'ZHONG ET AL. · 2025 ONLINE',title:'Honeybee agroecosystems',organism:'Apis mellifera · blueberry / cranberry',text:'20 matched tissue pools · proteins + RNA',n:20,views:2};
demoCatalog.splice(demoCatalog.findIndex(d=>d.id==='honeybee-eggs-2025'),1,beeEntry);
for(const tissue of ['abdomen','gut'])demoCatalog.push({...beeEntry,id:'bee-own-'+tissue+'-2025',path:'data/demos/bee-own-'+tissue+'-2025.json',hidden:true});
const renderExperienceBeforeBee=renderExperience;
renderExperience=function(){renderExperienceBeforeBee();const intro=document.querySelector('.demo-heading p');if(intro)intro.textContent='Use matched protein and RNA measurements from honeybee tissues, or compare the human datasets.';const q=new URLSearchParams(location.search),wanted=q.get('study');if(wanted==='bee'){const tissue=['head','abdomen','gut'].includes(q.get('tissue'))?q.get('tissue'):'head';loadDemo('bee-own-'+tissue+'-2025');}else if(wanted&&wanted!=='ccle'&&demoCatalog.some(d=>d.id===wanted))loadDemo(wanted);};
const loadDemoBeforeBee=loadDemo;
loadDemo=async function(id){
 const bee=id.startsWith('bee-own-');if(bee){pathwaySettings={view:'Proteomics',catalog:'GO_BP',mode:'ora',top_n:20,direction:'absolute'};if($('start-lab'))$('start-lab').disabled=true;}
 await loadDemoBeforeBee(id);
 if(activeDemo===id){const q=new URLSearchParams();q.set('study',bee?'bee':id);if(bee&&dataset.teaching.tissue!=='Head')q.set('tissue',id.includes('abdomen')?'abdomen':'gut');history.replaceState(null,'',location.pathname+'?'+q+'#study');}
 if(dataset?.teaching?.family==='bee-own-2025'){
  annotations={Proteomics:dataset.gene_labels||{},RNA:dataset.gene_labels||{}};
  document.querySelectorAll('[data-demo]').forEach(b=>{const own=b.dataset.demo===beeEntry.id;b.classList.toggle('selected',own);b.setAttribute('aria-pressed',String(own));b.querySelector('b').textContent=own?'Selected dataset':'Open dataset →';});
 }
 if($('start-lab')&&!busy)$('start-lab').disabled=false;
};
const summaryBeforeBee=renderStudySummary;
renderStudySummary=function(){
 summaryBeforeBee();$('bee-tissue-choice')?.remove();if(dataset?.teaching?.family!=='bee-own-2025')return;
 $('public-study-overview').insertAdjacentHTML('beforebegin',`<section id="bee-tissue-choice" class="bee-tissue-choice"><div><span class="kicker">ONE TISSUE. TWO MOLECULAR VIEWS.</span><h2>Choose a tissue</h2><p>Each tissue has 10 blueberry and 10 cranberry apiary-year pools. Switching tissue starts a new analysis.</p></div><div class="bee-tissue-buttons">${[['head','Head'],['abdomen','Abdomen'],['gut','Midgut']].map(([id,name])=>`<button data-bee-tissue="${id}" aria-pressed="${dataset.teaching.tissue===name}">${name}<small>20 paired pools</small></button>`).join('')}</div></section>`);
 document.querySelectorAll('[data-bee-tissue]').forEach(b=>b.onclick=()=>loadDemo('bee-own-'+b.dataset.beeTissue+'-2025'));
};
const busyBeforeBee=setBusy;
setBusy=function(value){busyBeforeBee(value);document.querySelectorAll('[data-bee-tissue]').forEach(b=>b.disabled=value);};
const pathwaysBeforeBee=renderPathwaysPanel;
renderPathwaysPanel=function(box){
 if(dataset.teaching?.annotation_mode!=='bee_genes'){pathwaysBeforeBee(box);return;}
 if(needMofa(box))return;
 const savedCatalog=pathwaySettings.catalog;pathwaysBeforeBee(box);
 box.querySelector('.explorer-banner p').textContent='Explore Apis mellifera GO biological processes and KEGG pathways using the matched bee gene IDs.';
 const update=()=>{$('path-catalog').innerHTML='<option value="GO_BP">Honeybee GO biological process</option><option value="Bee_KEGG">Honeybee KEGG pathways</option>';$('path-catalog').value=['GO_BP','Bee_KEGG'].includes(pathwaySettings.catalog)?pathwaySettings.catalog:'GO_BP';$('path-mode').value='ora';$('path-mode').querySelector('option[value="coverage"]').disabled=true;applySemanticControls();};
 pathwaySettings.catalog=['GO_BP','Bee_KEGG'].includes(savedCatalog)?savedCatalog:'GO_BP';update();$('path-view').onchange=update;
 box.insertAdjacentHTML('beforeend','<p class="notice">Mapping coverage is shown with each result. The background contains fitted bee genes with annotations. GO includes electronic annotations; KEGG uses Apis mellifera gene links.</p>');
};
const drawPathwaysBeforeBee=drawPathways;
drawPathways=function(e){drawPathwaysBeforeBee(e);if(dataset.teaching?.annotation_mode!=='bee_genes')return;const notes=document.querySelectorAll('#pathway-content > p.note');if(notes.length)notes[notes.length-1].textContent=pathwayMethodNote(e)+' Bee GO: UniProt taxon 7460, retrieved 2026-11-11, with available GO.db 3.19.1 ancestors. KEGG: official Apis mellifera (ame) gene pathways, retrieved 2026-11-11.';document.querySelector('#pathway-content > a[href="data/metabolite_mapping.csv"]')?.remove();};
const codeBeforeBee=currentRSnippet;
currentRSnippet=function(){
 if(dataset?.teaching?.family!=='bee-own-2025')return codeBeforeBee();
 if(current==='pathways')return `source("R/lab_explore.R")\nannotations <- jsonlite::fromJSON("bee_annotations.json", simplifyVector=FALSE)\nora <- factor_enrichment(lab, mofa, annotations,\n  factor=${Number($('path-factor')?.value??factorIndex)+1}, view="${$('path-view')?.value||'RNA'}",\n  catalog="${$('path-catalog')?.value||'GO_BP'}", top_n=${Number($('path-top')?.value??20)})\nplot_pathway_results(ora)`;
 if(['study','quality','single','integration','export'].includes(current))return `# In the extracted honeybee R kit:\n# Choose tissue <- "head", "abdomen" or "gut" in run_demo.R.\nsource("run_demo.R")\n# First installation: source("R/install_packages.R")\n# Try factors=2, 3, 4, 5 or 6.\n# Feature scaling: "unit_variance", "pareto", "none".\n# To reproduce source preparation: source("prepare_bee.R"); prepare_bee()`;
 return codeBeforeBee();
};
