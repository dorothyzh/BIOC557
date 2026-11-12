'use strict';
const courseGlossary=[

 ['Pareto','Pareto scaling','Divides each feature by the square root of its standard deviation; centering is controlled separately.'],
 ['ρ','Spearman rank correlation coefficient','Describes how two quantities vary together by rank.'],
 ['κ','Cohen kappa similarity','Measures agreement between term-gene membership profiles beyond chance.'],
 ['UniProt','Universal Protein Resource','A curated resource of protein sequences and functional annotations.'],
 ['MOFA','Multi-Omics Factor Analysis','Finds latent factors that describe variation across several molecular layers.'],
 ['MOFA2','Multi-Omics Factor Analysis, version 2','The R software used to fit and inspect MOFA models.'],
 ['PCA','Principal Component Analysis','Summarizes the main directions of variation in a matrix.'],
 ['PC','Principal Component','One axis obtained from principal component analysis.'],
 ['GO','Gene Ontology','A structured vocabulary describing gene functions and biological processes.'],
 ['BP','Biological Process','The Gene Ontology branch describing biological activities.'],
 ['GO_BP','Gene Ontology: Biological Process','The biological-process annotation catalog used here.'],
 ['KEGG','Kyoto Encyclopedia of Genes and Genomes','A database linking genes and compounds to pathways.'],
 ['ORA','Over-Representation Analysis','Tests whether an annotation occurs more often in a selected feature list than in its background.'],
 ['DAVID','Database for Annotation, Visualization and Integrated Discovery','A resource for gene annotation, enrichment and grouping related functional terms.'],
 ['EASE','Expression Analysis Systematic Explorer','The modified Fisher-test score used by DAVID.'],
 ['BH','Benjamini-Hochberg procedure','Adjusts p values for multiple comparisons.'],
 ['FDR','False Discovery Rate','The expected proportion of false discoveries among selected results under the testing procedure.'],
 ['q','BH-adjusted p value in this lab','A p value adjusted with the Benjamini-Hochberg procedure.'],
 ['p','p value','The probability of a result at least this extreme under the specified null model.'],
 ['R²','Coefficient of determination','Fraction of the model input reconstructed relative to a zero baseline; with centering on, it describes variation around the mean.'],
 ['K','Number of latent factors','The chosen number of factors in the model.'],
 ['SD','Standard Deviation','Measures the spread of a feature. Unit-variance scaling divides by SD; Pareto scaling divides by its square root.'],
 ['QC','Quality Control','Checks measurement gaps, scales and sample patterns.'],
 ['CPM','Counts Per Million','RNA counts divided by the effective library size, multiplied by one million.'],
 ['TMM','Trimmed Mean of M-values','An RNA library-size normalization method that trims extreme log-ratios.'],
 ['TPM','Transcripts Per Million','RNA abundance adjusted for transcript length and library composition.'],
 ['TMT','Tandem Mass Tag','An isotope-labeling method for multiplexed protein quantification.'],
 ['MS','Mass Spectrometry','Measures ions by mass-to-charge ratio.'],
 ['LC','Liquid Chromatography','Separates compounds before measurement.'],
 ['LC–MS','Liquid Chromatography-Mass Spectrometry','Combines compound separation with mass-spectrometric detection.'],
 ['MS/MS','Tandem Mass Spectrometry','Fragments ions to help identify molecules.'],
 ['m/z','Mass-to-charge ratio','The mass of an ion divided by its charge.'],
 ['RT','Retention Time','The time at which a compound elutes from a chromatography column.'],
 ['HBB','Highbush Blueberry','One crop environment in the honeybee study.'],
 ['CRA','Cranberry','The other crop environment in the honeybee study.'],
 ['CCLE','Cancer Cell Line Encyclopedia','A collection of molecular measurements from cancer cell-line models.'],
 ['T2D','Type 2 Diabetes','The disease group in the human-islet dataset.'],
 ['RRID','Research Resource Identifier','The source identifier used to match the islet donors.'],
 ['LOC','Locus identifier prefix','An NCBI label for a gene locus; its number is an identifier rather than a functional name.'],
 ['GeneID','NCBI Gene Identifier','A database identifier for a gene.'],
 ['WGCNA','Weighted Gene Co-expression Network Analysis','Groups genes with correlated expression patterns.'],
 ['HTSeq','High-Throughput Sequencing analysis framework','Software used to count RNA-sequencing reads.'],
 ['RNA-seq','RNA Sequencing','Measures RNA by sequencing RNA-derived molecules.'],
 ['ROS','Reactive Oxygen Species','Reactive oxygen-containing molecules involved in signaling and oxidative stress.'],
 ['TCA','Tricarboxylic Acid cycle','A central metabolic pathway that oxidizes carbon fuels.'],
 ['Dscam','Down syndrome cell adhesion molecule','A cell-adhesion gene with roles in neuronal wiring and extensive alternative splicing.'],
];
function glossaryMatches(text){return courseGlossary.filter(([s])=>{if(s==='LOC')return /\bLOC[0-9]+/.test(text);if(s==='GeneID')return /\bGeneID[_:]?[0-9]+/.test(text);if(s==='PC')return /\bPC(?:[0-9]+)?\b/.test(text);const escaped=s.replace(/[.*+?^${}()|[\]\\]/g,'\\$&');return new RegExp('(^|[^A-Za-z0-9])'+escaped+'([^A-Za-z0-9]|$)',['R','p','q','K','ID'].includes(s)?'':'i').test(text);});}
function glossaryMarkup(rows){return rows.map(([s,full,meaning])=>`<div><b>${esc(s)}</b><span>${esc(full)}. ${esc(meaning)}</span></div>`).join('');}
function updateCourseGlossary(){
 const hero=document.querySelector('.hero-copy>p');if(hero&&!document.querySelector('.hero-abbreviations'))hero.insertAdjacentHTML('afterend','<p class="hero-abbreviations">PCA: principal component analysis · MOFA: multi-omics factor analysis.</p>');
 // Clear every generated note before reading chart text, including old nested notes.
 document.querySelectorAll('.chart-abbreviations').forEach(note=>note.remove());
 const target=current==='study'?$('study-panel'):$('chapter-content');if(!target)return;
 const explained=new Set();
 target.querySelectorAll('.chart-section').forEach(section=>{
  if(section.querySelector('.chart-section')||!section.querySelector('.plot'))return;
  const rows=glossaryMatches(section.innerText).filter(([term])=>!explained.has(term));
  rows.forEach(([term])=>explained.add(term));
  if(rows.length)section.insertAdjacentHTML('beforeend',`<div class="chart-abbreviations">${rows.map(([s,f,m])=>`<span data-term="${esc(s)}"><b>${esc(s)}</b>: ${esc(f)} — ${esc(m)}</span>`).join('')}</div>`);
 });
 // Reuse the reference panel so redraws retain its open state and keyboard focus.
 let reference=$('course-glossary');
 if(!reference){reference=document.createElement('details');reference.id='course-glossary';reference.className='course-glossary';reference.innerHTML=`<summary>Key terms and abbreviations</summary><div>${glossaryMarkup(courseGlossary)}</div>`;}
 if(reference.parentElement!==target)target.append(reference);
}
const renderBeforeGlossary=renderChapter;renderChapter=function(){renderBeforeGlossary();updateCourseGlossary();};
const studyBeforeGlossary=renderStudySummary;renderStudySummary=function(){studyBeforeGlossary();updateCourseGlossary();};
let glossaryFrame=null;
function scheduleCourseGlossary(){if(glossaryFrame!==null)return;glossaryFrame=requestAnimationFrame(()=>{glossaryFrame=null;updateCourseGlossary();});}
const plotBeforeGlossary=plot;plot=function(...args){const p=plotBeforeGlossary(...args);Promise.resolve(p).then(scheduleCourseGlossary);return p;};

const setChapterBeforeGlossary=setChapter;setChapter=function(id){setChapterBeforeGlossary(id);updateCourseGlossary();};
