'use strict';
chapters.splice(1,0,['upload','Upload your data','02 / YOUR DATA','Bring your own molecular measurements.','Load at least two omics tables, check the sample IDs, then follow the same analysis steps.']);
chapters.forEach((chapter,index)=>{chapter[2]=String(index+1).padStart(2,'0')+' / '+chapter[2].split(' / ').slice(1).join(' / ');});
const experienceBeforeOwnData=renderExperience;
renderExperience=function(){experienceBeforeOwnData();if(!$('open-own-data')){const button=document.createElement('button');button.id='open-own-data';button.textContent='Upload your data →';button.onclick=()=>setChapter('upload');document.querySelector('.hero-actions').append(button);}};
const chapterBeforeOwnData=setChapter;
setChapter=function(id){chapterBeforeOwnData(id);$('show-code').hidden=id==='upload';};
let ownDataLoading=false;
const uploadViews=[['Proteomics','protein'],['RNA','rna'],['Metabolomics','metabolite']];
async function loadUpload(){
 if(busy||ownDataLoading)return;
 ownDataLoading=true;$('own-data-preview').hidden=true;$('load-upload').disabled=true;$('upload-message').textContent='Checking your files…';$('upload-message').className='upload-status';
 try{
  const selected=uploadViews.filter(([,key])=>$(key+'-file').files[0]),metadataFile=$('metadata-file').files[0];
  if(selected.length<2||!metadataFile)throw Error('Choose at least two omics tables and the metadata CSV.');
  const files=selected.map(([,key])=>$(key+'-file').files[0]);
  if([...files,metadataFile].reduce((sum,file)=>sum+file.size,0)>8*1024*1024)throw Error('Keep all selected files under 8 MB in total.');
  const parsed=await Promise.all([metadataFile,...files].map(loadCSV));
  const out=OwnDataFormat.study(parsed[0],selected.map(([name,key],i)=>({name,filename:files[i].name,parsed:parsed[i+1],scale:$(key+'-scale').value,orientation:$('upload-orientation').value})),{title:$('upload-title').value,organism:$('upload-organism').value});
  out.source.filenames.push(metadataFile.name);
  dataset=out;activeDemo='uploaded';annotations={};result=null;baseline=null;comparison=null;networkResult=null;enrichmentResult=null;runHistory.length=0;functionalSummary=null;summarySource=null;
  applyStudyDefaults(out);applyDemoFactorLimits();
  document.querySelectorAll('.rail-bottom a[download],#code-dialog a[download],#r-walkthrough-link').forEach(link=>link.hidden=true);
  document.querySelectorAll('[data-demo]').forEach(button=>{button.classList.remove('selected');button.setAttribute('aria-pressed','false');button.querySelector('b').textContent='Open dataset →';});
  $('public-study-overview').hidden=true;$('uploaded-overview').hidden=false;$('uploaded-title').textContent=out.title;
  $('uploaded-summary').textContent=out.sample_ids.length+' matched samples · '+out.views.map(v=>v.features.length+' '+v.name+' features').join(' · ');
  $('start-uploaded').onclick=()=>{setChapter('quality');runAnalysis();};
  const shared=out.views.find(v=>v.name==='Proteomics')?.features.filter(id=>out.views.find(v=>v.name==='RNA')?.features.includes(id))||[];
  $('own-data-preview').innerHTML=`<span class="kicker">READY FOR ANALYSIS</span><h2>${esc(out.title)}</h2><p>${out.sample_ids.length} matched samples · ${out.group_order.map(group=>esc(group)+' (n='+out.groups.filter(x=>x===group).length+')').join(' · ')}</p><div class="table-wrap"><table><thead><tr><th>Molecular view</th><th>Features</th><th>Missing values</th><th>File layout</th></tr></thead><tbody>${out.views.map(v=>`<tr><th>${esc(v.name)}</th><td>${v.features.length}</td><td>${pct(v.values.flat().filter(x=>x===null).length/(v.values.length*v.features.length))}</td><td>${esc(v.orientation)}</td></tr>`).join('')}</tbody></table></div><p class="note">${shared.length>=4?shared.length+' matching protein–RNA gene IDs found. The model aligns them in the same column order.':'Each view keeps its own feature IDs. Samples are aligned by sample_id.'} ${out.sample_ids.length<16?'Start with PCA for this small sample set.':'Start with MOFA + PCA using the selected preparation settings.'}</p><button class="primary" id="run-own-upload">Start analysis →</button>`;
  $('own-data-preview').hidden=false;$('run-own-upload').onclick=()=>{setChapter('quality');runAnalysis();};
  $('upload-message').textContent='Files checked. Samples are aligned and ready.';$('upload-message').className='upload-status success';
  history.replaceState(null,'',location.pathname+'#upload');setChapter('upload');
 }catch(error){$('upload-message').textContent=error.message;$('upload-message').className='upload-status error';}
 finally{ownDataLoading=false;$('load-upload').disabled=busy;}
};
const pathwayBeforeOwnData=renderPathwaysPanel;
renderPathwaysPanel=function(box){if(result.annotation_scope==='none'){box.innerHTML='<h2>Functional annotation</h2><p>This dataset uses the Other organism setting. Explore the factors, feature weights and network, or choose a supported gene catalog in <a href="#upload">Upload your data</a> and load the files again.</p>';return;}pathwayBeforeOwnData(box);};

for(const id of ['protein-file','rna-file','metabolite-file','metadata-file','protein-scale','rna-scale','metabolite-scale','upload-orientation','upload-organism','upload-title'])$(id).addEventListener('change',()=>{if(!$('own-data-preview').hidden){$('own-data-preview').hidden=true;$('upload-message').textContent='Files or settings changed. Select Check & load data to update the analysis input.';$('upload-message').className='upload-status';}});

document.querySelectorAll('[data-clear-file]').forEach(button=>button.onclick=event=>{event.preventDefault();event.stopPropagation();const input=$(button.dataset.clearFile);input.value='';input.dispatchEvent(new Event('change',{bubbles:true}));});

const demoBeforeOwnData=loadDemo;
loadDemo=async function(id){await demoBeforeOwnData(id);if(activeDemo!=='uploaded'){$('own-data-preview').hidden=true;$('upload-message').textContent='Select Check & load data to use your files.';$('upload-message').className='upload-status';}};
