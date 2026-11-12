let py=null, ready=null;
const send=(type,extra={})=>postMessage({type,...extra});
async function initialize(){
  if(ready)return ready;
  ready=(async()=>{
    send('progress',{message:'Starting your independent analysis workspace…'});
    importScripts('https://cdn.jsdelivr.net/pyodide/v0.29.3/full/pyodide.js');
    py=await loadPyodide({indexURL:'https://cdn.jsdelivr.net/pyodide/v0.29.3/full/',stdout:()=>{},stderr:()=>{}});
    send('progress',{message:'Loading scientific libraries. First visit takes longer; each learner has an independent workspace.'});
    await py.loadPackage(['numpy','pandas','scipy','scikit-learn','h5py','micropip']);
    const [code,wheel,exploreCode,catalogs,beeCatalogs,summaryCode,summaryRefs]=await Promise.all([fetch('lab_core.py?v=20261111g').then(r=>{if(!r.ok)throw Error('Analysis code could not load.');return r.text();}),fetch('vendor/mofapy2-0.7.5-py3-none-any.whl').then(r=>{if(!r.ok)throw Error('MOFA package could not load.');return r.arrayBuffer();}),fetch('lab_explore.py?v=20261111l').then(r=>r.text()),fetch('data/enrichment_annotations.json?v=20261111e').then(r=>r.text()),fetch('data/bee_enrichment_annotations.json?v=20261111g').then(r=>{if(!r.ok)throw Error('Bee annotations could not load.');return r.text();}),fetch('lab_summary.py?v=20261111l').then(r=>{if(!r.ok)throw Error('Summary code could not load.');return r.text();}),fetch('data/summary_references.json?v=20261111l').then(r=>{if(!r.ok)throw Error('Summary references could not load.');return r.text();})]);
    py.FS.writeFile('/home/pyodide/lab_summary.py',summaryCode);py.FS.writeFile('/home/pyodide/summary_references.json',summaryRefs);py.FS.writeFile('/home/pyodide/lab_core.py',code);py.FS.writeFile('/home/pyodide/lab_explore.py',exploreCode);py.FS.writeFile('/home/pyodide/enrichment_annotations.json',catalogs);py.FS.writeFile('/home/pyodide/bee_enrichment_annotations.json',beeCatalogs);
    py.FS.writeFile('/tmp/mofapy2-0.7.5-py3-none-any.whl',new Uint8Array(wheel));
    await py.runPythonAsync("import micropip\nawait micropip.install('emfs:/tmp/mofapy2-0.7.5-py3-none-any.whl', deps=False)\nfrom lab_core import run_json\nfrom lab_explore import compare_json, explore_json");
    send('ready');return py;
  })().catch(e=>{ready=null;send('error',{message:'The scientific environment could not start. '+e.message});throw e;});
  return ready;
}
onmessage=async({data})=>{
  if(data.type==='init'){await initialize().catch(()=>{});return;}
  if(data.type==='compare'||data.type==='explore'){
    try{await initialize();if(data.type==='compare'){py.globals.set('_dataset',JSON.stringify(data.dataset));py.globals.set('_parameters',JSON.stringify(data.parameters));py.globals.set('_counts',JSON.stringify(data.counts||[2,3,5]));}else{py.globals.set('_result',JSON.stringify(data.result));py.globals.set('_options',JSON.stringify(data.options||{}));}send('progress',{message:data.type==='compare'?'Fitting the selected factor counts…':'Calculating the selected network or pathway query…'});const value=await py.runPythonAsync(data.type==='compare'?'compare_json(_dataset, _parameters, _counts)':'explore_json(_result, _options)');send(data.type==='compare'?'comparison':'exploration',{requestId:data.requestId,result:JSON.parse(value),options:data.options});}catch(e){send('error',{requestId:data.requestId,message:String(e.message||e).split('\n').slice(-4).join('\n')});}return;
  }
  if(data.type==='run'){
    try{await initialize();send('progress',{message:data.parameters.method==='mofa'?'Fitting MOFA to your selected features…':'Calculating your analysis…'});py.globals.set('_dataset',JSON.stringify(data.dataset));py.globals.set('_parameters',JSON.stringify(data.parameters));const value=await py.runPythonAsync('run_json(_dataset, _parameters)');send('result',{requestId:data.requestId,result:JSON.parse(value)});}
    catch(e){send('error',{requestId:data.requestId,message:String(e.message||e).split('\n').slice(-4).join('\n')});}
  }
};
