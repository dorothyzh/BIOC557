# BIOC557: samples are rows; feature columns retain their original published scale.
# The same preparation decisions are used by the website. Students write only R.
read_lab_data <- function(root='.') {
  meta <- read.csv(file.path(root,'data/metadata.csv'), check.names=FALSE, stringsAsFactors=FALSE)
  files <- c(Proteomics='proteomics.csv', RNA='rna.csv', Metabolomics='metabolomics.csv')
  views <- lapply(files, function(f) {
    d <- read.csv(file.path(root,'data',f), check.names=FALSE, stringsAsFactors=FALSE)
    if (anyDuplicated(d$sample_id) || !setequal(d$sample_id,meta$sample_id)) stop('Sample IDs must match metadata exactly.')
    x <- as.matrix(d[match(meta$sample_id,d$sample_id),-1,drop=FALSE]); storage.mode(x)<-'double'
    rownames(x)<-meta$sample_id
    if (any(is.infinite(x)) || anyDuplicated(colnames(x))) stop('Invalid numeric values or duplicate features.')
    x
  })
  if(anyDuplicated(meta$sample_id) || length(unique(meta$group))!=2 || any(table(meta$group)<3)) stop('Use unique IDs and two groups with at least 3 samples each.')
  # Source table already contains log2 protein ratios, log2(TPM+1) RNA and log10 metabolites.
  list(metadata=meta, views=views, levels=unique(meta$group))
}
prepare_views <- function(study, max_missing=.30, model_features=100, standardize=TRUE, normalization='none',feature_scaling=NULL,mofa_center=TRUE) {
  if(is.null(feature_scaling))feature_scaling<-if(standardize)'unit_variance' else 'none'
  if(!feature_scaling%in%c('none','unit_variance','pareto'))stop('Unknown feature scaling.')
  lapply(study$views,function(raw) {
    x<-raw
    if(normalization=='sample_median') {m<-apply(x,1,median,na.rm=TRUE); if(any(!is.finite(m))) stop('A sample has no observed values.'); x<-x-m+median(m)}
    missing<-colMeans(is.na(x)); variance<-apply(x,2,var,na.rm=TRUE)
    keep<-missing<=max_missing & is.finite(variance) & variance>1e-12
    x<-x[,keep,drop=FALSE]; if(ncol(x)<4) stop('Too few variable features remain.')
    top<-order(-variance[keep], seq_len(sum(keep)))[seq_len(min(model_features,ncol(x)))]
    model<-x[,top,drop=FALSE]
    if(mofa_center)model<-sweep(model,2,colMeans(model,na.rm=TRUE))
    sd<-apply(model,2,stats::sd,na.rm=TRUE)
    if(feature_scaling=='unit_variance')model<-sweep(model,2,sd,FUN='/')
    if(feature_scaling=='pareto')model<-sweep(model,2,sqrt(sd),FUN='/')
    imputed<-model
    for(j in seq_len(ncol(imputed))) imputed[is.na(imputed[,j]),j]<-median(imputed[,j],na.rm=TRUE)
    list(raw=raw, retained=x, model=model, pca=imputed,
         input_features=ncol(raw), retained_features=ncol(x), missing_fraction=mean(is.na(raw)), imputed_entries=sum(is.na(model)))
  })
}
observed_welch <- function(x, groups, levels) {
  rows<-lapply(seq_len(ncol(x)),function(j){a<-x[groups==levels[1],j];b<-x[groups==levels[2],j];a<-a[is.finite(a)];b<-b[is.finite(b)]
    p<-if(length(a)>=3 && length(b)>=3) tryCatch(t.test(b,a,var.equal=FALSE)$p.value,error=function(e) NA_real_) else NA_real_
    data.frame(feature=colnames(x)[j],mean_difference=mean(b)-mean(a),p_value=p,observed_group_1=length(a),observed_group_2=length(b))})
  d<-do.call(rbind,rows);d$q_value<-NA_real_;ok<-is.finite(d$p_value);d$q_value[ok]<-p.adjust(d$p_value[ok],method='BH');d
}
stable_pca <- function(x) {
  p<-prcomp(x,center=TRUE,scale.=FALSE,rank.=3)
  for(j in seq_len(ncol(p$x))) {i<-which.max(abs(p$rotation[,j]));if(p$rotation[i,j]<0){p$x[,j]<- -p$x[,j];p$rotation[,j]<- -p$rotation[,j]}}
  list(scores=p$x,loadings=p$rotation,explained=p$sdev[1:3]^2/sum(p$sdev^2))
}
run_lab <- function(study,max_missing=.30,model_features=100,standardize=TRUE,normalization='none',feature_scaling=NULL,mofa_center=TRUE) {
  if(is.null(feature_scaling))feature_scaling<-if(standardize)'unit_variance' else 'none'
  prepared<-prepare_views(study,max_missing,model_features,standardize,normalization,feature_scaling,mofa_center)
  pcas<-lapply(prepared,function(v)stable_pca(v$pca))
  joint<-do.call(cbind,lapply(prepared,function(v)v$pca/sqrt(ncol(v$pca))))
  diff<-lapply(prepared,function(v)observed_welch(v$retained,study$metadata$group,study$levels))
  distances<-lapply(prepared,function(v)as.numeric(dist(v$pca)))
  agreement<-cor(do.call(cbind,distances),method='spearman')
  colnames(agreement)<-rownames(agreement)<-names(prepared)
  list(study=study,prepared=prepared,pca=pcas,joint_pca=stable_pca(joint),differential=diff,agreement=agreement,
       parameters=list(max_missing=max_missing,model_features=model_features,standardize=standardize,normalization=normalization,feature_scaling=feature_scaling,mofa_center=mofa_center))
}
fit_lab_mofa <- function(lab,factors=3,iterations=100,seed=557,runtime='managed',scale_views=FALSE,outfile=file.path('outputs','BIOC557_MOFA.hdf5')) {
  if(!requireNamespace('MOFA2',quietly=TRUE)) stop('Run R/install_packages.R first.')
  if(!runtime %in% c('managed','existing','basilisk')) stop('Choose a supported R calculation environment.')
  use_basilisk<-identical(runtime,'basilisk')
  if(runtime=='managed') {
    if(utils::packageVersion('reticulate') < '1.41') stop('Update reticulate in R, or choose runtime="basilisk".')
    Sys.setenv(RETICULATE_USE_MANAGED_VENV='yes')
    reticulate::py_require('mofapy2==0.7.5',python='3.12')
  }
  views<-lapply(names(lab$prepared),function(name){x<-t(lab$prepared[[name]]$model);rownames(x)<-paste(name,rownames(x),sep='::');x});names(views)<-names(lab$prepared)
  # Prefixes keep RNA and protein genes distinct inside the factor model.
  m<-MOFA2::create_mofa(views)
  data_options<-MOFA2::get_default_data_options(m);data_options$scale_views<-scale_views;data_options$scale_groups<-FALSE;data_options$center_groups<-if(is.null(lab$parameters$mofa_center))TRUE else lab$parameters$mofa_center;data_options$use_float32<-FALSE
  model_options<-MOFA2::get_default_model_options(m);model_options$num_factors<-factors;model_options$spikeslab_weights<-TRUE;model_options$ard_weights<-TRUE;model_options$likelihoods<-rep('gaussian',length(views))
  training_options<-MOFA2::get_default_training_options(m);training_options$maxiter<-iterations;training_options$seed<-seed;training_options$convergence_mode<-'fast';training_options$verbose<-FALSE;training_options$startELBO<-1;training_options$freqELBO<-5
  m<-MOFA2::prepare_mofa(m,data_options=data_options,model_options=model_options,training_options=training_options)
  dir.create(dirname(outfile),recursive=TRUE,showWarnings=FALSE)
  m<-MOFA2::run_mofa(m,outfile=outfile,use_basilisk=use_basilisk)
  z<-MOFA2::get_factors(m)[[1]];z<-z[lab$study$metadata$sample_id,,drop=FALSE]
  w<-MOFA2::get_weights(m)
  for(name in names(views)) w[[name]]<-w[[name]][rownames(views[[name]]),,drop=FALSE]
  for(j in seq_len(ncol(z))) {allw<-unlist(lapply(w,function(a)a[,j]));if(allw[which.max(abs(allw))]<0){z[,j]<- -z[,j];w<-lapply(w,function(a){a[,j]<- -a[,j];a})}}
  model_data<-MOFA2::get_data(m,add_intercept=FALSE,na.rm=FALSE)
  inputs<-lapply(names(views),function(name)t(model_data[[name]][[1]][rownames(views[[name]]),lab$study$metadata$sample_id,drop=FALSE]));names(inputs)<-names(views)
  r2<-do.call(rbind,lapply(names(views),function(name){x<-inputs[[name]];den<-sum(x^2,na.rm=TRUE);sapply(seq_len(ncol(z)),function(j)1-sum((x-outer(z[,j],w[[name]][,j]))^2,na.rm=TRUE)/den)}));rownames(r2)<-names(views);colnames(r2)<-colnames(z)
  total<-sapply(names(views),function(name){x<-inputs[[name]];1-sum((x-z%*%t(w[[name]]))^2,na.rm=TRUE)/sum(x^2,na.rm=TRUE)})
  list(model=m,scores=z,weights=w,variance_explained=r2,total_r2=total,model_input=inputs,center=data_options$center_groups,scale_views=scale_views)
}
