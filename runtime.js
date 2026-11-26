'use strict';
// Plot operations are tied to the actual DOM node, including during fast navigation.
const labGraphStates=new WeakMap(),labGraphNodes=new Set();
function releaseLabGraph(node){
 const state=labGraphStates.get(node);if(state)state.token++;
 labGraphNodes.delete(node);
 try{Plotly.purge(node);}catch{}
}
function renderLabPlot(id,traces,layout){
 const node=$(id);if(!node)return Promise.resolve(null);
 let state=labGraphStates.get(node);
 if(!state){state={token:0,pending:Promise.resolve()};labGraphStates.set(node,state);}
 const token=++state.token;labGraphNodes.add(node);
 state.pending=state.pending.catch(()=>null).then(async()=>{
  if(!node.isConnected||token!==state.token)return null;
  node.querySelector('.plot-error')?.remove();node.setAttribute('aria-busy','true');
  try{
   const graph=await (node.classList.contains('js-plotly-plot')?Plotly.react(node,traces,layout,plotConfig):Plotly.newPlot(node,traces,layout,plotConfig));
   if(!node.isConnected){releaseLabGraph(node);return null;}
   if(token!==state.token)return null;
   node.setAttribute('aria-busy','false');return graph;
  }catch(error){
   if(!node.isConnected||token!==state.token)return null;
   node.setAttribute('aria-busy','false');
   try{Plotly.purge(node);}catch{}
   const message=document.createElement('div');message.className='plot-error';message.setAttribute('role','alert');
   const text=document.createElement('p');text.textContent='This figure did not finish drawing. Your analysis is still available.';
   const retry=document.createElement('button');retry.textContent='Redraw figure';retry.onclick=()=>renderLabPlot(id,traces,layout);
   message.append(text,retry);node.replaceChildren(message);
   console.warn('Figure could not render:',id,error?.message||String(error));return null;
  }
 });
 return state.pending;
}
const labGraphObserver=new MutationObserver(()=>{
 for(const node of labGraphNodes)if(!node.isConnected)releaseLabGraph(node);
});
labGraphObserver.observe($('chapter-content'),{childList:true,subtree:true});
let labResizeFrame=null;
window.addEventListener('resize',()=>{
 if(labResizeFrame!==null)return;
 labResizeFrame=requestAnimationFrame(()=>{labResizeFrame=null;for(const node of labGraphNodes){if(node.isConnected&&node.getBoundingClientRect().width>0&&node.classList.contains('js-plotly-plot'))Promise.resolve(Plotly.Plots.resize(node)).catch(()=>{});}});
});

// Every network transfer has a deadline and one retry for temporary failures.
async function fetchLabAsset(url,options={}){
 let last;
 for(let attempt=0;attempt<2;attempt++){
  const controller=new AbortController(),timer=setTimeout(()=>controller.abort(),30000);
  try{const response=await fetch(url,{...options,signal:controller.signal});if(!response.ok)throw Error('Could not load '+url.split('?')[0]+' (HTTP '+response.status+').');
   // Consume the body while the deadline still applies.
   const bytes=await response.arrayBuffer();return new Response(bytes,{status:response.status,headers:response.headers});
  }catch(error){last=error;}finally{clearTimeout(timer);}
 }
 throw Error(last?.name==='AbortError'?'Loading took too long. Select Retry to try again.':last?.message||'The file could not load.');
}
let labEngineEpoch=0,labInitTimer=null,labRequestTimer=null;
function clearLabTimers(){clearTimeout(labInitTimer);clearTimeout(labRequestTimer);labInitTimer=labRequestTimer=null;}
function failLabRequest(message,restart=false){
 clearLabTimers();requestId++;
 if(restart){labEngineEpoch++;worker?.terminate();worker=null;ready=false;}
 setBusy(false);$('run-status').hidden=true;$('runtime-label').textContent=ready?'Analysis ready':'Ready to retry';
 $('runtime-dot').classList.toggle('ready',ready);$('error-banner').textContent=message;$('error-banner').hidden=false;$('retry-engine').hidden=false;
}
function armLabRequest(){
 clearTimeout(labRequestTimer);
 labRequestTimer=setTimeout(()=>failLabRequest('This run took longer than expected. Select Run analysis to restart it; the last completed results are still available.',true),240000);
}
initWorker=function(){
 clearLabTimers();worker?.terminate();ready=false;
 const epoch=++labEngineEpoch;
 $('retry-engine').hidden=true;$('runtime-label').textContent='Loading analysis engine';$('runtime-dot').classList.remove('ready');
 const active=new Worker('engine.worker.js?v=20261126a');worker=active;
 labInitTimer=setTimeout(()=>{if(epoch===labEngineEpoch)failLabRequest('The analysis engine could not finish loading. Select Retry to load it again.',true);},120000);
 active.onmessage=({data})=>{
  if(epoch!==labEngineEpoch)return;
  if(data.requestId!==undefined&&data.requestId!==requestId)return;
  if(data.type==='progress'){$('runtime-label').textContent=ready?'Calculating analysis':'Loading analysis engine';if(busy){$('run-status').textContent=data.message;$('run-status').hidden=false;}return;}
  if(data.type==='ready'){clearTimeout(labInitTimer);ready=true;$('runtime-label').textContent=busy?'Calculating analysis':'Analysis ready';$('runtime-dot').classList.add('ready');return;}
  if(data.type==='error'){failLabRequest(data.message,!ready);return;}
  if(['result','comparison','exploration'].includes(data.type)){
   clearTimeout(labRequestTimer);
   if(data.type==='result')finishRun(data.result);
   if(data.type==='comparison')finishComparison(data.result);
   if(data.type==='exploration')finishExploration(data.result,data.options);
  }
 };
 active.onerror=event=>{if(epoch===labEngineEpoch)failLabRequest(event.message||'The analysis engine stopped. Select Retry to restart it.',true);};
 active.postMessage({type:'init'});
};
const busyBeforeReliability=setBusy;
setBusy=function(value){busyBeforeReliability(value);$('cancel-analysis').hidden=!value;document.querySelectorAll('[data-demo],[data-bee-tissue],#run-empty').forEach(node=>node.disabled=value||labStudyLoading);};
$('cancel-analysis').onclick=()=>failLabRequest('Run stopped. Select Run analysis when you are ready to continue.',true);
async function retryLabLoading(){
 if(busy)failLabRequest('Previous run stopped. Select Run analysis to try again.',true);
 $('error-banner').hidden=true;
 if(!dataset){
  $('retry-engine').disabled=true;
  try{await boot();}catch(error){$('error-banner').textContent=error.message;$('error-banner').hidden=false;$('retry-engine').hidden=false;}
  finally{$('retry-engine').disabled=false;}
 }else initWorker();
}
$('retry-engine').onclick=retryLabLoading;
const runBeforeReliability=runAnalysis;
runAnalysis=async function(){if(busy||labStudyLoading)return;if(!dataset){$('error-banner').textContent='Wait for the selected dataset to finish loading.';$('error-banner').hidden=false;return;}if(!worker&&!useServer)initWorker();armLabRequest();return runBeforeReliability();};
const requestBeforeReliability=scientificRequest;
scientificRequest=async function(type,extra={}){if(busy)return;if(!worker&&!useServer)initWorker();armLabRequest();return requestBeforeReliability(type,extra);};
$('run-empty').onclick=()=>runAnalysis();

if(location.hostname.endsWith('.chatgpt.site')){
 const notice=document.createElement('p');notice.className='mirror-notice';
 notice.append('Course link: ');const link=document.createElement('a');link.href='https://dorothyzh.github.io/BIOC557/';link.textContent='dorothyzh.github.io/BIOC557';notice.append(link);
 document.querySelector('main').prepend(notice);
}
