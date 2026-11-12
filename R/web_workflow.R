# Functions used by the complete, student-facing chapter scripts.
web_install_packages <- function(use_mofa=TRUE) {
  old <- getOption('BIOC557.install_mofa')
  on.exit(options(BIOC557.install_mofa=old), add=TRUE)
  options(BIOC557.install_mofa=use_mofa)
  source('R/install_packages.R', local=TRUE)
  source('R/lab_core.R', local=environment(web_install_packages))
  source('R/lab_plots.R', local=environment(web_install_packages))
  source('R/lab_explore.R', local=environment(web_install_packages))
}
web_metadata <- function(input) {
  ids <- unlist(input$sample_ids, use.names=FALSE)
  if (is.null(input$metadata)) return(data.frame(sample_id=ids, group=unlist(input$groups)))
  keys <- unique(unlist(lapply(input$metadata, names)))
  rows <- lapply(input$metadata, function(row) {
    row <- setNames(lapply(keys, function(key) if(is.null(row[[key]])) NA else row[[key]]), keys)
    as.data.frame(row, stringsAsFactors=FALSE, check.names=FALSE)
  })
  meta <- do.call(rbind, rows)
  stopifnot(!anyDuplicated(meta$sample_id), setequal(meta$sample_id, ids))
  meta[match(ids, meta$sample_id), , drop=FALSE]
}
web_numeric_matrix <- function(view) {
  x <- do.call(rbind, lapply(view$values, function(row) vapply(row,
    function(value) if(is.null(value)) NA_real_ else as.numeric(value), numeric(1))))
  dimnames(x) <- list(unlist(view$sample_ids), unlist(view$features))
  stopifnot(!anyDuplicated(rownames(x)), !anyDuplicated(colnames(x)), !any(is.infinite(x)))
  x
}
web_prepare_input <- function(input, parameters, plan=list(kind='none')) {
  meta <- web_metadata(input); ids <- meta$sample_id
  views <- setNames(lapply(input$views, web_numeric_matrix),
                    vapply(input$views, function(view) view$name, character(1)))
  if(parameters$challenge=='shuffle_rows') {
    stopifnot(!is.null(plan$shuffle_operations))
    for(name in names(views)) {
      op <- plan$shuffle_operations[[name]]
      views[[name]] <- views[[name]][unlist(op$rows)+1L, unlist(op$columns)+1L, drop=FALSE]
    }
  }
  views <- lapply(views, function(x) {
    stopifnot(setequal(rownames(x), ids))
    x[match(ids, rownames(x)), , drop=FALSE]
  })
  if(parameters$challenge=='break_pairing') {
    stopifnot(length(plan$permutation)==length(ids))
    last <- length(views)
    views[[last]] <- views[[last]][unlist(plan$permutation)+1L, , drop=FALSE]
    rownames(views[[last]]) <- ids
  }
  if(parameters$challenge=='missing') {
    stopifnot(!is.null(plan$masked_indices_by_view))
    for(name in names(views)) {
      indices <- plan$masked_indices_by_view[[name]]
      if(length(indices)) views[[name]][do.call(rbind,lapply(indices,unlist))+1L] <- NA_real_
    }
  }
  units <- setNames(lapply(input$views, function(view) view$unit), names(views))
  for(view in input$views) {
    name <- view$name
    log_input <- parameters$transform=='log2' ||
                 (parameters$transform=='auto' && view$input_scale=='linear')
    if(log_input) {
      stopifnot(view$input_scale=='linear', all(views[[name]][is.finite(views[[name]])]>0))
      views[[name]] <- log2(views[[name]])
      units[[name]] <- 'log2 abundance'
    }
  }
  list(metadata=meta, views=views, units=units,
       levels=if(length(input$group_order)) unlist(input$group_order) else unique(meta$group))
}
web_save_plot <- function(figure, name, width=10, height=6) {
  figure <- figure + ggplot2::theme(panel.grid=ggplot2::element_blank(),
    plot.title=ggplot2::element_text(face='bold',color='#111820'))
  print(figure)
  ggplot2::ggsave(paste0('outputs/',name,'.png'),figure,width=width,height=height,dpi=300,bg='white')
  ggplot2::ggsave(paste0('outputs/',name,'.pdf'),figure,width=width,height=height,bg='white')
  invisible(figure)
}
web_pca_plot <- function(lab, view=NULL) {
  p <- if(is.null(view))lab$joint_pca else lab$pca[[view]]
  scatter_scores(p$scores,lab,if(is.null(view))'Joint PCA' else paste(view,'PCA'),
    sprintf('PC1 (%.1f%%)',100*p$explained[1]),sprintf('PC2 (%.1f%%)',100*p$explained[2]))
}
web_quality_plots <- function(lab, view) {
  x <- lab$prepared[[view]]$raw; values <- lab$prepared[[view]]$retained
  samples <- rownames(x); groups <- lab$study$metadata$group
  missing <- data.frame(sample=factor(samples,levels=samples),missing=rowMeans(is.na(x)),group=groups)
  box <- data.frame(sample=factor(rep(samples,ncol(values)),levels=samples),value=as.numeric(values),group=rep(groups,ncol(values)))
  colors <- setNames(c('#287D9E','#BD683C'),lab$study$levels)
  list(missing=ggplot2::ggplot(missing,ggplot2::aes(sample,missing,fill=group))+
    ggplot2::geom_col()+ggplot2::scale_fill_manual(values=colors)+
    ggplot2::scale_y_continuous(labels=scales::percent,limits=c(0,1))+
    ggplot2::labs(title='Missing measurements',x='Sample',y='Fraction missing')+lab_theme()+
    ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55,hjust=1)),
    distribution=ggplot2::ggplot(box,ggplot2::aes(sample,value,fill=group))+
      ggplot2::geom_boxplot(outlier.shape=NA,na.rm=TRUE)+ggplot2::scale_fill_manual(values=colors)+
      ggplot2::labs(title='Abundance after preparation',x='Sample',y=lab$study$units[[view]])+lab_theme()+
      ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55,hjust=1)),pca=web_pca_plot(lab,view))
}
web_volcano <- function(lab, view) {
  d <- lab$differential[[view]]; d$significant <- is.finite(d$q_value)&d$q_value<.05
  ggplot2::ggplot(d,ggplot2::aes(mean_difference,-log10(pmax(q_value,1e-300)),color=significant))+
    ggplot2::geom_point(size=1.7,alpha=.8,na.rm=TRUE)+
    ggplot2::geom_hline(yintercept=-log10(.05),linetype=2,color='#b9a078')+
    ggplot2::scale_color_manual(values=c('FALSE'='#96b7b9','TRUE'='#bd683c'),guide='none')+
    ggplot2::labs(title=paste(view,'group comparison'),x=paste(lab$study$levels[2],'minus',lab$study$levels[1]),y='-log10(BH q)')+lab_theme()
}
web_weights <- function(lab, mofa, factor=1, n=8) {
  d <- do.call(rbind,lapply(names(mofa$weights),function(view) {
    w <- mofa$weights[[view]][,factor]; ix <- head(order(-abs(w)),n)
    data.frame(view=view,feature=colnames(lab$prepared[[view]]$model)[ix],weight=w[ix])
  }))
  d$key <- paste(d$view,d$feature,sep='::'); labels <- setNames(d$feature,d$key)
  d$key <- factor(d$key,levels=d$key[order(d$view,d$weight)])
  ggplot2::ggplot(d,ggplot2::aes(weight,key,fill=view))+ggplot2::geom_col()+
    ggplot2::facet_wrap(~view,scales='free_y',ncol=1)+ggplot2::scale_y_discrete(labels=labels)+
    ggplot2::scale_fill_manual(values=view_colors,guide='none')+
    ggplot2::labs(title=paste('Factor',factor,'features'),x='MOFA feature weight',y='Feature')+lab_theme()
}
web_associations <- function(lab,mofa) {
  meta <- lab$study$metadata; traits <- list()
  for(key in intersect(c('group','age','sex','year','province','rna_library_size','rna_TMM_factor'),names(meta))) {
    value <- meta[[key]]; ok <- !is.na(value)&value!=''
    if(sum(ok)<4 || length(unique(value[ok]))<2)next
    numeric <- suppressWarnings(as.numeric(as.character(value)))
    if(all(is.finite(numeric[ok])))traits[[key]] <- numeric
    else if(length(unique(value[ok]))==2) {
      levels <- unique(value[ok]);traits[[paste0(key,': ',levels[2],' (1)')]] <- match(value,levels)-1L
    }
  }
  for(view in names(lab$prepared)) {
    value <- rowMeans(is.finite(mofa$model_input[[view]]))
    if(length(unique(value))>1)traits[[paste(view,'observed fraction')]] <- value
  }
  rows <- lapply(names(traits),function(name)do.call(rbind,lapply(seq_len(ncol(mofa$scores)),function(k) {
    x <- traits[[name]];y <- mofa$scores[,k];ok <- is.finite(x)&is.finite(y)
    data.frame(trait=name,factor=paste('Factor',k),rho=if(sum(ok)>=4)cor(x[ok],y[ok],method='spearman') else NA_real_,n=sum(ok))
  })))
  do.call(rbind,rows)
}
web_association_plot <- function(d) {
  ggplot2::ggplot(d,ggplot2::aes(factor,trait,fill=rho))+ggplot2::geom_tile()+
    ggplot2::geom_text(ggplot2::aes(label=ifelse(is.finite(rho),sprintf('%.2f',rho),'')),size=3)+
    ggplot2::scale_fill_gradient2(low='#287d9e',mid='white',high='#bd683c',limits=c(-1,1))+
    ggplot2::labs(title='Factors and sample characteristics',x='MOFA factor',y='Sample characteristic',fill='Spearman rho')+lab_theme()
}
web_feature_correlations <- function(lab,mofa,factor=1) {
  values <- do.call(cbind,lapply(names(mofa$weights),function(view) {
    ix <- head(order(-abs(mofa$weights[[view]][,factor])),4)
    x <- mofa$model_input[[view]][,ix,drop=FALSE]
    colnames(x) <- paste(view,colnames(lab$prepared[[view]]$model)[ix],sep=': ');x
  }))
  rho <- suppressWarnings(cor(values,use='pairwise.complete.obs',method='spearman'))
  observed <- crossprod(is.finite(values));rho[observed<4] <- NA_real_
  d <- as.data.frame(as.table(rho));names(d) <- c('first','second','rho')
  ggplot2::ggplot(d,ggplot2::aes(first,second,fill=rho))+ggplot2::geom_tile()+
    ggplot2::scale_fill_gradient2(low='#287d9e',mid='white',high='#bd683c',limits=c(-1,1))+
    ggplot2::labs(title='Selected features across layers',x='Molecular feature',y='Molecular feature',fill='Spearman rho')+
    lab_theme()+ggplot2::theme(axis.text.x=ggplot2::element_text(angle=55,hjust=1,size=8),axis.text.y=ggplot2::element_text(size=8))
}
web_pca_model <- function(lab) {
  p <- lab$joint_pca; offset <- 0L; weights <- inputs <- list();r2 <- list()
  for(view in names(lab$prepared)) {
    x <- lab$prepared[[view]]$pca/sqrt(ncol(lab$prepared[[view]]$pca))
    x <- sweep(x,2,colMeans(x)); n <- ncol(x)
    w <- p$loadings[offset+seq_len(n),,drop=FALSE];offset <- offset+n
    rownames(w) <- colnames(x);weights[[view]] <- w;inputs[[view]] <- x
    r2[[view]] <- sapply(seq_len(ncol(p$scores)),function(k)1-sum((x-outer(p$scores[,k],w[,k]))^2)/sum(x^2))
  }
  list(scores=p$scores,weights=weights,model_input=inputs,variance_explained=do.call(rbind,r2))
}
web_group_terms <- function(rows, universe, annotations, category) {
  if(!nrow(rows))return(data.frame())
  sets <- lapply(seq_len(nrow(rows)),function(i) {
    terms <- annotations$catalogs[[rows$catalog[i]]]
    term <- terms[[which(vapply(terms,function(t)t$id,character(1))==rows$id[i])]]
    unlist(term$members,use.names=FALSE)
  })
  distance <- matrix(0,nrow(rows),nrow(rows));N <- length(universe)
  for(i in seq_len(nrow(rows)))if(i>1)for(j in seq_len(i-1)) {
    similarity <- 0
    if(length(intersect(sets[[i]],sets[[j]]))>=3 && N>0) {
      a <- intersect(sets[[i]],universe);b <- intersect(sets[[j]],universe)
      po <- (length(intersect(a,b))+length(setdiff(universe,union(a,b))))/N
      pe <- (length(a)*length(b)+(N-length(a))*(N-length(b)))/N^2
      if(pe<1)similarity <- (po-pe)/(1-pe)
    }
    distance[i,j] <- distance[j,i] <- 1-max(0,similarity)
  }
  active <- seq_len(nrow(rows));members <- lapply(active,identity)
  while(length(active)>1) {
    pairs <- which(upper.tri(distance),arr.ind=TRUE)
    pairs <- pairs[pairs[,1]%in%active & pairs[,2]%in%active,,drop=FALSE]
    values <- round(distance[pairs],12)
    choice <- order(values,pairs[,1],pairs[,2])[1]
    if(values[choice]>.5)break
    i <- pairs[choice,1];j <- pairs[choice,2]
    for(other in setdiff(active,c(i,j)))distance[i,other] <- distance[other,i] <- max(distance[i,other],distance[j,other])
    members[[i]] <- c(members[[i]],members[[j]]);active <- setdiff(active,j)
  }
  groups <- integer(nrow(rows))
  for(i in active)groups[members[[i]]] <- i
  do.call(rbind,lapply(unique(groups),function(group) {
    d <- rows[groups==group,,drop=FALSE]
    d <- d[order(d$report_q,d$p_value,d$term_in_background,d$name),,drop=FALSE]
    genes <- unique(unlist(strsplit(d$features,';',fixed=TRUE)))
    data.frame(factor=d$factor[1],category=category(d$name[1]),leading_term=d$name[1],
      report_q=min(d$report_q),score=mean(-log10(pmax(d$p_value,1e-300))),
      terms=nrow(d),selected_features=length(genes),views=paste(sort(unique(d$view)),collapse=' + '))
  }))
}
web_functional_summary <- function(lab,model,annotations,top_n=20,categories=NULL,annotation_mode='') {
  queries <- list(); all <- list(); kmax <- ncol(model$scores)
  for(k in seq_len(kmax))for(view in names(lab$prepared)) {
    if(model$variance_explained[view,k]<.01)next
    catalogs <- if(view=='Metabolomics') {
      if(annotation_mode=='human_genes_only')character() else 'Metabolite_pathways'
    } else intersect(c('GO_BP','Bee_KEGG','Reactome'),names(annotations$catalogs))
    for(catalog in catalogs)for(direction in c('positive','negative')) {
      q <- factor_enrichment(lab,model,annotations,k,view,catalog,top_n,direction)
      if(!nrow(q$rows))next
      q$rows$factor <- k;q$rows$view <- view;q$rows$catalog <- catalog;q$rows$direction <- direction
      all[[length(all)+1L]] <- q$rows
    }
  }
  tests <- if(length(all))do.call(rbind,all) else data.frame()
  if(nrow(tests))tests$report_q <- p.adjust(tests$p_value,'BH')
  get_category <- function(name) {
    if(length(categories))for(item in categories)if(grepl(item$pattern,name,ignore.case=TRUE,perl=TRUE))return(item$title)
    if(grepl('metabolic|catabolic|biosynth',name,ignore.case=TRUE))'Metabolism and biosynthesis' else 'Functional annotation profile'
  }
  themes <- list()
  factors <- lapply(seq_len(kmax),function(k) {
    rows <- if(nrow(tests))tests[tests$factor==k,,drop=FALSE] else data.frame()
    pool <- if(nrow(rows))rows[rows$overlap>=2 & rows$overlap*rows$background_size>rows$query_size*rows$term_in_background & rows$term_in_background<=.6*rows$background_size,,drop=FALSE] else rows
    if(nrow(pool))pool <- pool[order(pool$report_q,pool$p_value,-pool$overlap),,drop=FALSE]
    candidates <- head(pool[pool$p_value<=.1,,drop=FALSE],80)
    if(!nrow(candidates))candidates <- head(pool,8)
    gene_universe <- unique(unlist(lapply(lab$prepared[intersect(c('Proteomics','RNA'),names(lab$prepared))],function(v)colnames(v$model))))
    compound_universe <- if('Metabolomics'%in%names(lab$prepared))colnames(lab$prepared$Metabolomics$model) else character()
    grouped <- rbind(web_group_terms(candidates[candidates$view!='Metabolomics',,drop=FALSE],gene_universe,annotations,get_category),
      web_group_terms(candidates[candidates$view=='Metabolomics',,drop=FALSE],compound_universe,annotations,get_category))
    if(nrow(grouped))grouped <- head(grouped[order(grouped$report_q,-grouped$score),,drop=FALSE],5)
    themes[[k]] <<- grouped
    leading <- if(nrow(grouped))grouped[1,,drop=FALSE] else NULL
    data.frame(factor=k,category=if(is.null(leading))'Functional assignment open' else leading$category,
      leading_term=if(is.null(leading))'' else leading$leading_term,
      report_q=if(is.null(leading))NA_real_ else leading$report_q,
      significant_tests=if(nrow(rows))sum(rows$report_q<.05) else 0,
      views_above_1pct=sum(model$variance_explained[,k]>=.01))
  })
  list(factors=do.call(rbind,factors),tests=tests,themes=do.call(rbind,themes),contributions=model$variance_explained)
}
web_summary_plot <- function(summary,component='Factor') {
  f <- summary$factors; kmax <- nrow(f);d <- as.data.frame(as.table(summary$contributions))
  names(d) <- c('view','factor','r2');d$factor <- rep(seq_len(kmax),each=nrow(summary$contributions))
  views <- rownames(summary$contributions);d$y <- kmax-d$factor+1+(match(d$view,views)-(length(views)+1)/2)*.17
  f$y <- kmax-f$factor+1; f$label <- paste(component,f$factor)
  f$status <- ifelse(is.finite(f$report_q),paste0(ifelse(f$report_q<.05,'BH-supported','Candidate'),' | q = ',formatC(f$report_q,format='g',digits=3)),'No grouped annotation')
  ggplot2::ggplot()+ggplot2::geom_segment(data=d,ggplot2::aes(x=0,xend=pmax(0,r2),y=y,yend=y),color='#e7efee',linewidth=5)+
    ggplot2::geom_segment(data=d,ggplot2::aes(x=0,xend=pmax(0,r2),y=y,yend=y,color=view),linewidth=4)+
    ggplot2::geom_text(data=d,ggplot2::aes(x=-.025,y=y,label=view),hjust=1,size=3)+
    ggplot2::geom_text(data=d,ggplot2::aes(x=pmax(0,r2)+.025,y=y,label=sprintf('%.1f%%',100*r2)),hjust=0,size=3)+
    ggplot2::geom_point(data=f,ggplot2::aes(x=.95,y=y),size=13,color='#193f4e')+
    ggplot2::geom_text(data=f,ggplot2::aes(x=.95,y=y,label=paste0(if(component=='PC')'PC' else 'F',factor)),color='white',size=3.5)+
    ggplot2::geom_text(data=f,ggplot2::aes(x=1.08,y=y+.08,label=category),hjust=0,fontface='bold',size=3.5)+
    ggplot2::geom_text(data=f,ggplot2::aes(x=1.08,y=y-.15,label=status),hjust=0,size=3)+
    ggplot2::scale_color_manual(values=view_colors)+ggplot2::coord_cartesian(xlim=c(-.2,2.3),ylim=c(.45,kmax+.6),clip='off')+
    ggplot2::labs(title='Factor functions and evidence',subtitle='View contributions and report-wide BH-adjusted enrichment',color='Molecular view')+
    ggplot2::theme_void()+ggplot2::theme(legend.position='bottom',plot.title=ggplot2::element_text(face='bold',size=16,color='#111820'),
      plot.subtitle=ggplot2::element_text(size=10),plot.margin=ggplot2::margin(18,25,18,25))
}
