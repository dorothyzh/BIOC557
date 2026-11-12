'use strict';
function downloadSummaryPdf(s,save=true){
 if(!s||typeof jspdf==='undefined')return;
 const {jsPDF}=jspdf,doc=new jsPDF({unit:'mm',format:'a4',compress:true});const W=210,H=297,M=18,B=277;let y=22,X=M,C=W-2*M;
 const ink=[28,58,70],teal=[31,116,128],muted=[88,109,116];
 const ascii=x=>String(x??'').replace(/<[^>]*>/g,'').replace(/[\u2010-\u2015\u2212]/g,'-').replace(/\u00a0/g,' ').replace(/[‘’]/g,"'").replace(/[“”]/g,'"').replace(/ρ/g,'rho').replace(/κ/g,'kappa').replace(/²/g,'2').replace(/₁/g,'1').replace(/₂/g,'2').replace(/₀/g,'0').replace(/≥/g,'>=').replace(/≤/g,'<=').replace(/↔/g,' / ').replace(/→/g,' to ').replace(/[^\x20-\x7e\n\u00a1-\u00ff]/g,' ');
 const line=(x1,y1,x2,y2)=>{doc.setDrawColor(210,225,224);doc.line(x1,y1,x2,y2);};
 function newPage(){doc.addPage();y=24;}
 function ensure(h){if(y+h>B)newPage();}
 function text(str,{size=10,bold=false,color=ink,gap=3,width=C}={}){doc.setFont('helvetica',bold?'bold':'normal');doc.setFontSize(size);doc.setTextColor(...color);const lines=doc.splitTextToSize(ascii(str),width),lh=size*.45;for(const t of lines){ensure(lh);doc.text(t,X,y);y+=lh;}y+=gap;}
 function heading(str){ensure(16);text(str,{size:16,bold:true,color:teal,gap:4});}
 function small(str){text(str,{size:8.5,color:muted,gap:3});}
 function refLine(r){text(`${r.title}. ${r.journal}, ${r.year}.`,{size:9,gap:1});if(r.url){ensure(6);doc.setFontSize(8);doc.setTextColor(...teal);doc.textWithLink(ascii(r.url),M,y,{url:r.url});y+=7;}}
 doc.setProperties({title:'BIOC557 Functional Summary - '+s.dataset,subject:'Current analysis, functional results and supporting literature',author:'BIOC557 | Huan Zhong'});
 small('BIOC557 / MULTI-OMICS LAB');text('Functional summary',{size:26,bold:true,gap:6});text(s.dataset,{size:16,bold:true,gap:5});
 small(`${s.model} | ${s.samples} matched samples | ${s.views.join(' + ')} | Workspace ${sessionId}`);
 text(`Prepared from the completed browser run. ${s.factors.length} components; top ${s.top_n} features per sign and view.`,{size:11});
 const p=s.parameters;small(`Feature alignment: ${p.feature_alignment}. Feature scaling: ${p.feature_scaling}. Center: ${p.mofa_center?'on':'off'}. Whole-view scaling: ${p.mofa_scale_views?'on':'off'}. Maximum missingness: ${Math.round(p.max_missing*100)}%. Seed: ${p.seed}.`);
 heading('Component overview');
 for(const f of s.factors){ensure(22);text(`${f.label} | ${f.themes[0]?.category||'Functional assignment open'}`,{size:11,bold:true,gap:1});small(f.contributions.map(v=>`${v.view}: ${(100*v.r2).toFixed(1)}% reconstruction`).join(' | '));small(`${f.total_significant_terms} term tests pass report-wide BH q < 0.05. ${f.themes[0]?.status||'Open assignment'}.`);}
 small('Reconstruction basis: '+s.reconstruction_basis);
 text('Reading the report',{size:11,bold:true});small('BH-supported themes contain a term with report-wide BH q < 0.05. Candidate themes are ranked functional leads. Each factor has at most one interpretation supported by enrichment and reviewed primary studies. Positive and negative weights describe opposite ends of a component; the sign can flip between fitted models.');
 const allRefs=[...s.study_papers,...s.references].filter((r,i,a)=>a.findIndex(x=>x.url===r.url)===i),refIds=Object.fromEntries(s.references.map(r=>[r.id,allRefs.findIndex(x=>x.url===r.url)+1]));
 for(const f of s.factors){newPage();X=M;C=W-2*M;small('COMPONENT PROFILE');heading(`${f.label}: ${f.themes[0]?.category||'Functional assignment open'}`);
  for(const v of f.contributions){ensure(9);doc.setFontSize(9);doc.setTextColor(...ink);doc.text(ascii(v.view),M,y);doc.setFillColor(231,239,238);doc.rect(M+34,y-3,105,4,'F');doc.setFillColor(...teal);doc.rect(M+34,y-3,105*Math.max(0,Math.min(1,v.r2)),4,'F');doc.text(`${(100*v.r2).toFixed(1)}%`,M+144,y);y+=7;}
  y+=4;const topY=y;C=82;
  text('Functional themes',{size:11,bold:true});
  if(!f.themes.length)small('Functional assignment is open. Examine the leading measured features and sample associations.');
  for(const g of f.themes.slice(0,3)){text(`${g.category} | ${g.status}`,{size:9,bold:true,gap:1});small(`${g.title}. Score ${g.score.toFixed(2)}; ${g.terms.length} related terms; lowest report BH q ${g.report_q.toPrecision(2)}.`);
   for(const t of g.terms.slice(0,2))text(`${t.name} [${t.id}]: ${t.view}, ${t.direction}; hits ${t.overlap}/${t.query_size}; report q=${t.report_q.toPrecision(2)}.`,{size:8,gap:2});
   text('IDs: '+g.genes.slice(0,5).join(', '),{size:8,color:muted,gap:4});
  }
  if(f.associations.length){text('Sample associations',{size:10,bold:true});small(f.associations.slice(0,3).map(a=>`${a.trait}: rho=${a.rho.toFixed(2)}, n=${a.n}`).join('; '));}
  const leftEnd=y;X=110;y=topY;C=82;text(f.hypotheses.length?'Main interpretation and evidence':'Supporting features',{size:11,bold:true});
  for(const h of f.hypotheses){text(h.title,{size:9,bold:true,gap:1});text(h.text,{size:9,gap:2});small('Data: '+h.basis);for(const p of h.paper_support)small('['+refIds[p.id]+'] '+p.support);text('Validation: '+h.test,{size:8.5,gap:2});y+=2;}
  text('Leading measured features',{size:10,bold:true});for(const [view,rows] of Object.entries(f.strongest_features))small(`${view}: `+rows.slice(0,4).map(r=>`${r.label} (${r.weight.toFixed(3)})`).join(', '));
  X=M;C=W-2*M;y=Math.max(y,leftEnd);
 }
 newPage();heading('Methods and reproducibility');text(s.method,{size:10});small(`Report-wide correction covers ${s.tested_term_queries} term tests across components, views, annotation catalogs and weight directions. Complete-linkage term clustering uses kappa >= ${s.cluster_kappa}, with at least three shared model genes. Positive/negative queries select up to ${s.top_n} features each. Views contributing at least 1% reconstruction enter annotation queries.`);
 text('Interpretation selection',{size:12,bold:true});small(s.interpretation_policy);text('Settings',{size:12,bold:true});for(const [key,value] of Object.entries(s.parameters))small(`${key}: ${value}`);
 newPage();heading('References');for(const [i,r] of allRefs.entries()){ensure(23);refLine({...r,title:'['+(i+1)+'] '+r.title});}
 newPage();heading('Abbreviations and explanations');const glossaryTop=y;X=M;C=82;let column=0;
 for(const [term,full,explanation] of courseGlossary){const lines=doc.splitTextToSize(ascii(full+' '+explanation),82).length;const need=10+lines*3.5;if(y+need>B){if(column===0){column=1;X=110;y=glossaryTop;}else{newPage();X=M;heading('Abbreviations, continued');column=0;}}text(`${term} - ${full}`,{size:9,bold:true,gap:1});small(explanation);}
 X=M;C=W-2*M;
 const pages=doc.getNumberOfPages();for(let i=1;i<=pages;i++){doc.setPage(i);line(M,285,W-M,285);doc.setFont('helvetica','normal');doc.setFontSize(8);doc.setTextColor(...muted);doc.text('BIOC557 | '+ascii(s.model)+' functional summary',M,291);doc.text(`${i} / ${pages}`,W-M,291,{align:'right'});}
 const bytes=doc.output('arraybuffer');if(save){doc.save('BIOC557_Functional_Summary_'+sessionId+'.pdf');if($('summary-status'))$('summary-status').textContent='PDF downloaded with the current settings, component profiles, evidence and references.';}return bytes;
}
