'use strict';
const OwnDataFormat=(()=>{
 const clean=value=>String(value??'').replace(/^\uFEFF/,'').trim();
 const preview=items=>items.slice(0,4).join(', ')+(items.length>4?'…':'');
 function table(parsed,label){
  if(Object.keys(parsed.meta?.renamedHeaders||{}).length)throw Error(label+': duplicate column names. Give every column a unique name.');
  const raw=parsed.meta?.fields||[],headers=raw.map(clean);
  if(headers.length<2||headers.some(x=>!x)||new Set(headers).size!==headers.length)throw Error(label+': use at least two columns with unique, nonempty names.');
  return{headers,rows:parsed.data.map(row=>Object.fromEntries(raw.map((key,i)=>[headers[i],clean(row[key])])))};
 }
 function metadata(parsed){
  const t=table(parsed,'Metadata');
  if(!t.headers.includes('sample_id')||!t.headers.includes('group'))throw Error('Metadata needs columns named sample_id and group.');
  if(t.rows.length<6||t.rows.length>300)throw Error('Use 6–300 samples, with at least 3 samples in each of two groups.');
  const ids=t.rows.map(row=>row.sample_id);
  if(ids.some(id=>!id))throw Error('Metadata has an empty sample_id.');
  if(new Set(ids).size!==ids.length)throw Error('Metadata has duplicate sample IDs. Each sample needs one row.');
  const groups=t.rows.map(row=>row.group),levels=[...new Set(groups)];
  if(groups.some(x=>!x)||levels.length!==2||levels.some(g=>groups.filter(x=>x===g).length<3))throw Error('Use exactly two groups, with at least 3 samples in each group.');
  return{ids,groups,levels,rows:t.rows};
 }
 function abundance(parsed,meta,name,scale,orientation='auto'){
  const t=table(parsed,name),idColumn=t.headers[0],sampleSet=new Set(meta.ids),columns=t.headers.slice(1);
  const columnSamples=columns.length===meta.ids.length&&columns.every(x=>sampleSet.has(x));
  const rowSamples=t.rows.length===meta.ids.length&&t.rows.every(row=>sampleSet.has(row[idColumn]));
  const transposed=orientation==='features_rows'||(orientation==='auto'&&columnSamples&&!rowSamples);
  let ids,features,raw;
  if(transposed){ids=columns;features=t.rows.map(row=>row[idColumn]);raw=ids.map(id=>t.rows.map(row=>row[id]));}
  else{ids=t.rows.map(row=>row[idColumn]);features=columns;raw=t.rows.map(row=>features.map(feature=>row[feature]));}
  if(new Set(ids).size!==ids.length)throw Error(name+': duplicate sample IDs.');
  const missing=meta.ids.filter(id=>!ids.includes(id)),extra=ids.filter(id=>!sampleSet.has(id));
  if(missing.length||extra.length)throw Error(name+': sample IDs must match metadata.'+(missing.length?' Missing: '+preview(missing)+'.':'')+(extra.length?' Extra: '+preview(extra)+'.':''));
  if(features.length<4||features.length>5000)throw Error(name+': use 4–5,000 feature IDs.');
  if(features.some(x=>!x)||new Set(features).size!==features.length)throw Error(name+': feature IDs must be unique and nonempty.');
  const values=raw.map((row,i)=>row.map((value,j)=>{
   if(value===''||/^(NA|N\/A|NaN)$/i.test(value))return null;
   const number=Number(value);
   if(!Number.isFinite(number))throw Error(name+': nonnumeric value at sample '+ids[i]+', feature '+features[j]+'. Use a number or a blank cell.');
   if(scale==='linear'&&number<=0)throw Error(name+': positive linear abundance must be greater than zero. Use blank cells for missing values, or choose Already transformed for log values.');
   return number;
  }));
  const order=meta.ids.map(id=>ids.indexOf(id)),aligned=order.map(i=>values[i]);
  if(aligned.some(row=>row.every(x=>x===null)))throw Error(name+': a sample has no observed measurements.');
  const variable=features.filter((_,j)=>new Set(aligned.map(row=>row[j]).filter(x=>x!==null)).size>1);
  if(variable.length<4)throw Error(name+': at least 4 features must vary across samples.');
  return{name,sample_ids:[...meta.ids],features,values:aligned,input_scale:scale,unit:scale==='linear'?'Uploaded positive abundance':'Uploaded transformed abundance',orientation:transposed?'Features in rows':'Samples in rows'};
 }
 function study(parsedMetadata,inputs,{title='Your multi-omics data',organism='none'}={}){
  if(inputs.length<2)throw Error('Choose at least two omics tables and the metadata CSV.');
  const meta=metadata(parsedMetadata),views=inputs.map(input=>abundance(input.parsed,meta,input.name,input.scale,input.orientation));
  const species={human:'Homo sapiens',bee:'Apis mellifera',none:'User-supplied organism'};
  return{title:clean(title)||'Your multi-omics data',sample_ids:meta.ids,groups:meta.groups,group_order:meta.levels,metadata:meta.rows,views,
   teaching:{id:'uploaded',species:species[organism]||species.none,annotation_mode:organism==='bee'?'bee_genes':organism==='human'?'human_genes':'none'},
   source:{kind:'user-upload',filenames:inputs.map(x=>x.filename||x.name),annotation_choice:organism}};
 }
 return{study,metadata,abundance};
})();
