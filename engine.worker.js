let py=null,ready=null,mofaReady=null,summaryReady=null;
const catalogsReady=new Map();
const send=(type,extra={})=>postMessage({type,...extra});
async function asset(url,binary=false){
 let last;
 for(let attempt=0;attempt<2;attempt++){
  const controller=new AbortController(),timer=setTimeout(()=>controller.abort(),30000);
  try{const response=await fetch(url,{signal:controller.signal});if(!response.ok)throw Error('HTTP '+response.status+' loading '+url.split('?')[0]);return binary?await response.arrayBuffer():await response.text();}
  catch(error){last=error;}finally{clearTimeout(timer);}
 }
 throw Error(last?.name==='AbortError'?'A required analysis file took too long to load. Select Retry.':last?.message||'An analysis file could not load.');
}
async function initialize(){
 if(ready)return ready;
 ready=(async()=>{
  send('progress',{message:'Loading the browser analysis engine…'});
  importScripts('https://cdn.jsdelivr.net/pyodide/v0.29.3/full/pyodide.js');
  py=await loadPyodide({indexURL:'https://cdn.jsdelivr.net/pyodide/v0.29.3/full/',stdout:()=>{},stderr:()=>{}});
  send('progress',{message:'Loading numerical libraries. First visit takes longer.'});
  const [code,exploreCode]=await Promise.all([asset('lab_core.py?v=20261126a'),asset('lab_explore.py?v=20261126a'),py.loadPackage(['numpy','scipy'])]);
  py.FS.writeFile('/home/pyodide/lab_core.py',code);py.FS.writeFile('/home/pyodide/lab_explore.py',exploreCode);
  await py.runPythonAsync('from lab_core import run_json\nfrom lab_explore import compare_json, explore_json');
  send('ready');return py;
 })().catch(error=>{ready=null;throw error;});
 return ready;
}
async function ensureMofa(){
 if(mofaReady)return mofaReady;
 mofaReady=(async()=>{
  send('progress',{message:'Loading MOFA libraries for this factor-model run…'});
  const [wheel]=await Promise.all([asset('vendor/mofapy2-0.7.5-py3-none-any.whl',true),py.loadPackage(['pandas','scikit-learn','h5py','micropip'])]);
  py.FS.writeFile('/tmp/mofapy2-0.7.5-py3-none-any.whl',new Uint8Array(wheel));
  await py.runPythonAsync("import micropip\nawait micropip.install('emfs:/tmp/mofapy2-0.7.5-py3-none-any.whl', deps=False)");
 })().catch(error=>{mofaReady=null;throw error;});
 return mofaReady;
}
async function ensureCatalogs(scope){
 if(scope==='none')return;
 const name=scope==='bee'?'bee_enrichment_annotations.json':'enrichment_annotations.json';
 if(catalogsReady.has(name))return catalogsReady.get(name);
 const pending=asset('data/'+name+'?v=20261111g').then(text=>py.FS.writeFile('/home/pyodide/'+name,text)).catch(error=>{catalogsReady.delete(name);throw error;});
 catalogsReady.set(name,pending);return pending;
}
async function ensureSummary(){
 if(summaryReady)return summaryReady;
 summaryReady=Promise.all([asset('lab_summary.py?v=20261126a'),asset('data/summary_references.json?v=20261126a')]).then(([code,refs])=>{py.FS.writeFile('/home/pyodide/lab_summary.py',code);py.FS.writeFile('/home/pyodide/summary_references.json',refs);}).catch(error=>{summaryReady=null;throw error;});
 return summaryReady;
}
async function handleMessage(data){
 try{
  await initialize();if(data.type==='init')return;
  if(data.type==='run'){
   if(data.parameters.method==='mofa')await ensureMofa();
   send('progress',{requestId:data.requestId,message:data.parameters.method==='mofa'?'Fitting MOFA to your selected features…':'Calculating PCA and group comparisons…'});
   py.globals.set('_dataset',JSON.stringify(data.dataset));py.globals.set('_parameters',JSON.stringify(data.parameters));
   const value=await py.runPythonAsync('run_json(_dataset, _parameters)');send('result',{requestId:data.requestId,result:JSON.parse(value)});
  }else if(data.type==='compare'){
   await ensureMofa();send('progress',{requestId:data.requestId,message:'Fitting the selected factor counts…'});
   py.globals.set('_dataset',JSON.stringify(data.dataset));py.globals.set('_parameters',JSON.stringify(data.parameters));py.globals.set('_counts',JSON.stringify(data.counts||[2,3,5]));
   const value=await py.runPythonAsync('compare_json(_dataset, _parameters, _counts)');send('comparison',{requestId:data.requestId,result:JSON.parse(value)});
  }else if(data.type==='explore'){
   const action=data.options?.action;
   if(action==='enrichment'||action==='summary'){
    send('progress',{requestId:data.requestId,message:'Loading '+(data.result.annotation_scope==='bee'?'honeybee':'human')+' functional annotations…'});
    await ensureCatalogs(data.result.annotation_scope);
   }
   if(action==='summary')await ensureSummary();
   send('progress',{requestId:data.requestId,message:action==='summary'?'Calculating the functional summary…':'Calculating the selected network or pathway query…'});
   py.globals.set('_result',JSON.stringify(data.result));py.globals.set('_options',JSON.stringify(data.options||{}));
   const value=await py.runPythonAsync('explore_json(_result, _options)');send('exploration',{requestId:data.requestId,result:JSON.parse(value),options:data.options});
  }
 }catch(error){send('error',{...(data.type==='init'?{}:{requestId:data.requestId}),message:String(error.message||error).split('\n').slice(-4).join('\n')});}
}
// Pyodide has a single interpreter. Never execute two Python jobs concurrently.
let taskQueue=Promise.resolve();
onmessage=({data})=>{taskQueue=taskQueue.then(()=>handleMessage(data));};
