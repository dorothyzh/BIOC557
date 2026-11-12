# Helpers for the short, beginner-facing default report.
source('R/lab_core.R')
source('R/lab_plots.R')
source('R/lab_explore.R')
read_course_study <- function() {
  config<-jsonlite::fromJSON('lab_config.json',simplifyVector=FALSE)
  data<-jsonlite::fromJSON(config$study_file,simplifyVector=FALSE)
  metadata<-do.call(rbind,lapply(data$metadata,as.data.frame,stringsAsFactors=FALSE))
  stopifnot(!anyDuplicated(metadata$sample_id))
  views<-setNames(lapply(data$views,function(view){
    x<-do.call(rbind,lapply(view$values,function(row)vapply(row,function(v)if(is.null(v))NA_real_ else as.numeric(v),numeric(1))))
    colnames(x)<-unlist(view$features);rownames(x)<-unlist(view$sample_ids)
    stopifnot(!anyDuplicated(rownames(x)),!anyDuplicated(colnames(x)),setequal(rownames(x),metadata$sample_id))
    x[match(metadata$sample_id,rownames(x)),,drop=FALSE]
  }),vapply(data$views,function(v)v$name,character(1)))
  settings<-jsonlite::fromJSON('default_settings.json',simplifyVector=FALSE)[[config$id]]
  list(config=config,data=data,settings=settings,study=list(metadata=metadata,views=views,levels=unlist(data$group_order),units=setNames(lapply(data$views,function(v)v$unit),names(views))))
}
course_settings_table <- function(course) {
  p<-course$settings
  data.frame(Setting=c('Dataset','Matched samples','Analysis','MOFA factors','Model features per view','Feature scaling','Protein/RNA alignment','Center features','Scale whole views','Seed','Maximum iterations'),
    Value=c(course$config$title,nrow(course$study$metadata),if(p$method=='mofa')'Multi-Omics Factor Analysis (MOFA)' else 'Principal Component Analysis (PCA)',if(p$method=='mofa')p$factors else 'Not used',p$model_features,p$feature_scaling,'Same gene IDs; same column order','Yes','No','557',if(p$method=='mofa')p$iterations else 'Not used'))
}
course_qc_table <- function(lab) do.call(rbind,lapply(names(lab$prepared),function(name){v<-lab$prepared[[name]];data.frame(View=name,Samples=nrow(v$raw),Input_features=ncol(v$raw),Retained_features=ncol(v$retained),Model_features=ncol(v$model),Missing_percent=round(100*mean(is.na(v$raw)),1))}))
course_save_figure <- function(figure,name,width=9,height=5.8) {
  dir.create('outputs/figures',recursive=TRUE,showWarnings=FALSE)
  ggplot2::ggsave(file.path('outputs/figures',paste0(name,'.png')),figure,width=width,height=height,dpi=300,bg='white')
  ggplot2::ggsave(file.path('outputs/figures',paste0(name,'.pdf')),figure,width=width,height=height,device=grDevices::pdf,useDingbats=FALSE,bg='white')
  ggplot2::ggsave(file.path('outputs/figures',paste0(name,'.svg')),figure,width=width,height=height,device=svglite::svglite,bg='white')
  figure
}
course_pca_panel <- function(lab) {
  plots<-lapply(names(lab$pca),function(view){p<-lab$pca[[view]];scatter_scores(p$scores,lab,paste(view,'PCA'),sprintf('PC1 (%.1f%%)',100*p$explained[1]),sprintf('PC2 (%.1f%%)',100*p$explained[2]))})
  patchwork::wrap_plots(plots,ncol=if(length(plots)>2)2 else 2,guides='collect')
}
course_joint_pca <- function(lab) {
  p<-lab$joint_pca;scatter_scores(p$scores,lab,'Joint PCA',sprintf('PC1 (%.1f%%)',100*p$explained[1]),sprintf('PC2 (%.1f%%)',100*p$explained[2]))
}
course_weight_plot <- function(lab,mofa,factor=1) {
  d<-do.call(rbind,lapply(names(mofa$weights),function(view){w<-mofa$weights[[view]][,factor];ix<-head(order(-abs(w)),8);ids<-colnames(lab$prepared[[view]]$model)[ix];data.frame(view=view,feature=ids,weight=as.numeric(w[ix]),key=paste(view,ids,sep='::'))}))
  d$key<-factor(d$key,levels=d$key[order(d$view,d$weight)]);labels<-setNames(d$feature,as.character(d$key))
  ggplot2::ggplot(d,ggplot2::aes(weight,key,fill=view))+ggplot2::geom_col(width=.7)+ggplot2::facet_wrap(~view,scales='free_y',ncol=1)+ggplot2::scale_y_discrete(labels=labels)+ggplot2::scale_fill_manual(values=view_colors,guide='none')+ggplot2::labs(title=paste('Factor',factor,': strongest features in each view'),x='MOFA feature weight',y=NULL)+lab_theme()
}
course_save_results <- function(course,lab,mofa=NULL) {
  dir.create('outputs',showWarnings=FALSE)
  for(view in names(lab$prepared)) {
    x<-lab$prepared[[view]]$model
    write.csv(data.frame(sample_id=rownames(x),x,check.names=FALSE),file.path('outputs',paste0(view,'_model_matrix.csv')),row.names=FALSE)
    write.csv(lab$differential[[view]],file.path('outputs',paste0(view,'_differential.csv')),row.names=FALSE)
  }
  genes<-colnames(lab$prepared$Proteomics$model)
  stopifnot(identical(genes,colnames(lab$prepared$RNA$model)))
  write.csv(data.frame(gene_id=genes,protein_column=genes,rna_column=genes),'outputs/gene_alignment.csv',row.names=FALSE)
  write.csv(course$study$metadata,'outputs/sample_metadata.csv',row.names=FALSE)
  write.csv(lab$joint_pca$scores,'outputs/joint_pca_scores.csv')
  if(!is.null(mofa)) {write.csv(mofa$scores,'outputs/mofa_scores.csv');write.csv(mofa$variance_explained,'outputs/mofa_factor_r2.csv')}
  saveRDS(list(lab=lab,mofa=mofa,config=course$config,settings=course$settings),'outputs/default_analysis.rds')
  record<-list(dataset=course$config$id,settings=course$settings,seed=557,center=TRUE,scale_views=FALSE,matched_sample_ids=course$study$metadata$sample_id,aligned_gene_ids=genes,protein_rna_columns_identical=TRUE,input_md5=unname(tools::md5sum(course$config$study_file)),software=lapply(c('R'=as.character(getRversion()),rmarkdown=as.character(packageVersion('rmarkdown')),ggplot2=as.character(packageVersion('ggplot2'))),identity))
  if(!is.null(mofa))record$software$MOFA2<-as.character(packageVersion('MOFA2'))
  jsonlite::write_json(record,'outputs/run_manifest.json',auto_unbox=TRUE,pretty=TRUE)
  writeLines(capture.output(sessionInfo()),'outputs/sessionInfo.txt')
  invisible(record)
}
course_references <- function(course) {
  papers<-course$data$source$papers
  if(is.null(papers))return(invisible(NULL))
  lines<-vapply(papers,function(p)paste0('- [',p$title,'](',p$url,'). *',p$journal,'*, ',p$year,'.'),character(1))
  knitr::asis_output(paste(lines,collapse='\n\n'))
}
