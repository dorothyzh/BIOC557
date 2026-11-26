# Integrated multi-omics graph: observed correlations, matching gene IDs and functions.
network_layout <- function(nodes,edges,view_names) {
  N <- nrow(nodes);anchors <- positions <- matrix(0,N,2);functions <- nodes$kind=='function'
  for(i in seq_along(view_names)) {
    angle <- pi*.8-(i-1)*2*pi/max(2,length(view_names));center <- c(1.45*cos(angle),1.25*sin(angle)-.35)
    indices <- which(nodes$view==view_names[i])
    for(k in seq_along(indices)) {
      j <- indices[k];theta <- k*2.399963229728653+(i-1)*.4;radius <- .22+.42*sqrt((k-.5)/length(indices))
      anchors[j,] <- center;positions[j,] <- center+radius*c(cos(theta),sin(theta))
    }
  }
  terms <- which(functions)
  for(k in seq_along(terms)) {
    columns <- min(4,length(terms));row <- (k-1)%/%columns;col <- (k-1)%%columns
    anchors[terms[k],] <- positions[terms[k],] <- c((col-(columns-1)/2)*1.02,1.65+row*.66)
  }
  for(iteration in seq_len(180)) {
    dx <- outer(positions[,1],positions[,1],'-');dy <- outer(positions[,2],positions[,2],'-');den <- (dx^2+dy^2+.04)^1.5
    force <- .004*cbind(rowSums(dx/den),rowSums(dy/den))+(anchors-positions)*ifelse(functions,.24,.055)
    for(e in seq_len(nrow(edges))) {
      i <- match(edges$source[e],nodes$id);j <- match(edges$target[e],nodes$id);delta <- positions[j,]-positions[i,];len <- max(1e-8,sqrt(sum(delta^2)))
      strength <- if(edges$kind[e]=='annotation').005 else .008;spring <- strength*(len-.55)*delta/len
      force[i,] <- force[i,]+spring;force[j,] <- force[j,]-spring
    }
    positions <- positions+pmax(-.12,pmin(.12,force))*.75
  }
  nodes$x <- positions[,1];nodes$y <- positions[,2];nodes
}
integrated_network <- function(lab,model,annotations=list(catalogs=list()),factor=NULL,top_per_view=20,
                               threshold=.6,q_max=.05,term_limit=8,include_within=FALSE,include_identity=TRUE) {
  if(!is.null(factor)&&(!factor%in%seq_len(ncol(model$scores))))stop('Choose an existing component.')
  if(top_per_view<4||top_per_view>40||term_limit<0||term_limit>12)stop('Choose 4-40 features and 0-12 function nodes.')
  nodes <- list();values <- list()
  for(view in names(lab$prepared)) {
    w <- model$weights[[view]];r2 <- pmax(0,model$variance_explained[view,])
    if(is.null(factor)) {
      weighted <- sweep(sweep(abs(w),2,pmax(1e-12,apply(abs(w),2,max)),'/'),2,sqrt(r2),'*')
      importance <- rowSums(weighted);strongest <- max.col(weighted,ties.method='first')
    } else {importance <- abs(w[,factor]);strongest <- rep(factor,nrow(w))}
    scale <- max(1e-12,max(importance));indices <- head(order(-importance,seq_along(importance)),top_per_view)
    for(j in indices) {
      k <- strongest[j];id <- colnames(lab$prepared[[view]]$model)[j]
      nodes[[length(nodes)+1L]] <- data.frame(id=paste(view,id,sep='::'),kind='feature',feature=id,label=id,view=view,
        column=j,weight=w[j,k],factor=k,factor_r2=r2[k],importance=importance[j]/scale,catalog='',q_value=NA_real_,p_value=NA_real_,status='',memberships='',views='',source_url='',overlap=NA_integer_,stringsAsFactors=FALSE)
      values[[length(values)+1L]] <- lab$prepared[[view]]$model[,j]
    }
  }
  nodes <- do.call(rbind,nodes);tests <- list()
  for(i in seq_len(nrow(nodes)-1L))for(j in (i+1L):nrow(nodes)) {
    if(!include_within&&nodes$view[i]==nodes$view[j])next
    x <- values[[i]];y <- values[[j]];ok <- is.finite(x)&is.finite(y);n <- sum(ok)
    if(n<10||sd(x[ok])<1e-12||sd(y[ok])<1e-12)next
    rho <- cor(x[ok],y[ok],method='spearman')
    p <- 2*pt(-abs(rho)*sqrt((n-2)/max(1e-30,1-rho^2)),df=n-2)
    tests[[length(tests)+1L]] <- data.frame(source=nodes$id[i],target=nodes$id[j],value=rho,p_value=p,n=n,kind='correlation',view_pair=paste(sort(c(nodes$view[i],nodes$view[j])),collapse=' / '),q_value=NA_real_,gene='',catalog='')
  }
  empty_edges <- data.frame(source=character(),target=character(),value=double(),p_value=double(),n=integer(),kind=character(),view_pair=character(),q_value=double(),gene=character(),catalog=character())
  correlation_tests <- if(length(tests))do.call(rbind,tests) else empty_edges
  correlation_tests$q_value <- p.adjust(correlation_tests$p_value,'BH')
  eligible <- correlation_tests[abs(correlation_tests$value)>=threshold & correlation_tests$q_value<=q_max,,drop=FALSE]
  eligible <- eligible[order(-abs(eligible$value),eligible$q_value,eligible$source,eligible$target,method='radix'),,drop=FALSE]
  chosen <- integer();pairs <- unique(eligible$view_pair)
  for(pair in pairs)chosen <- c(chosen,head(which(eligible$view_pair==pair),max(1,160%/%max(1,length(pairs)))))
  chosen <- head(c(chosen,setdiff(seq_len(nrow(eligible)),chosen)),160)
  edges <- eligible[chosen,,drop=FALSE]
  if(include_identity&&all(c('Proteomics','RNA')%in%names(lab$prepared))) {
    genes <- intersect(colnames(lab$prepared$Proteomics$model),colnames(lab$prepared$RNA$model))
    for(gene in genes) {
      a <- paste('Proteomics',gene,sep='::');b <- paste('RNA',gene,sep='::')
      if(all(c(a,b)%in%nodes$id))edges <- rbind(edges,data.frame(source=a,target=b,value=1,p_value=NA_real_,n=NA_integer_,kind='matched_gene',view_pair='',q_value=NA_real_,gene=gene,catalog=''))
    }
  }
  term_tests <- list()
  if(term_limit>0)for(view in names(lab$prepared)) {
    catalogs <- if(view=='Metabolomics')intersect('Metabolite_pathways',names(annotations$catalogs)) else intersect(c('GO_BP','Reactome','Bee_KEGG'),names(annotations$catalogs))
    chosen_ids <- nodes$feature[nodes$view==view]
    for(catalog in catalogs) {
      terms <- annotations$catalogs[[catalog]];annotated <- unique(unlist(lapply(terms,function(t)unlist(t$members,use.names=FALSE))))
      universe <- intersect(colnames(lab$prepared[[view]]$model),annotated);query <- intersect(chosen_ids,universe);N <- length(universe);n <- length(query)
      for(term in terms) {
        members <- intersect(unlist(term$members,use.names=FALSE),universe);M <- length(members)
        if(M<if(view=='Metabolomics')2 else 3)next
        hits <- sort(intersect(query,members));k <- length(hits)
        term_tests[[length(term_tests)+1L]] <- data.frame(id=term$id,name=term$name,catalog=catalog,view=view,p_value=if(n)phyper(k-1,M,N-M,n,lower.tail=FALSE) else 1,overlap=k,background_size=N,query_size=n,term_in_background=M,features=paste(hits,collapse=';'),source_url=if(is.null(term$source_url))'' else term$source_url)
      }
    }
  }
  functional_tests <- if(length(term_tests))do.call(rbind,term_tests) else data.frame()
  if(nrow(functional_tests))functional_tests$q_value <- p.adjust(functional_tests$p_value,'BH')
  candidates <- list()
  if(nrow(functional_tests)) {
    pool <- functional_tests[functional_tests$overlap>=2 & (functional_tests$catalog=='Metabolite_pathways' | functional_tests$term_in_background<=.6*functional_tests$background_size),,drop=FALSE]
    for(key in unique(paste(pool$catalog,pool$id,sep='::'))) {
      rows <- pool[paste(pool$catalog,pool$id,sep='::')==key,,drop=FALSE]
      best <- rows[order(rows$q_value,rows$p_value,rows$term_in_background,rows$name,method='radix')[1],,drop=FALSE]
      term <- annotations$catalogs[[best$catalog]][[which(vapply(annotations$catalogs[[best$catalog]],function(t)t$id==best$id,logical(1)))[1]]]
      memberships <- sort(nodes$id[nodes$feature %in% unlist(term$members,use.names=FALSE) & ((nodes$view=='Metabolomics')==(best$catalog=='Metabolite_pathways'))])
      candidates[[length(candidates)+1L]] <- list(id=paste0('term::',key),label=best$name,feature=best$id,catalog=best$catalog,q_value=min(rows$q_value),p_value=min(rows$p_value),views=sort(unique(nodes$view[nodes$id %in% memberships])),memberships=memberships,source_url=best$source_url,queries=rows)
    }
  }
  if(length(candidates))candidates <- candidates[order(vapply(candidates,function(t)t$q_value,numeric(1)),vapply(candidates,function(t)t$p_value,numeric(1)),-vapply(candidates,function(t)length(t$views),integer(1)),-vapply(candidates,function(t)length(t$memberships),integer(1)),vapply(candidates,function(t)t$label,character(1)),method='radix')]
  picked <- list()
  add <- function(candidate) {
    if(length(picked)>=term_limit||any(vapply(picked,function(t)t$id==candidate$id,logical(1))))return()
    for(old in picked)if(old$catalog==candidate$catalog&&length(intersect(old$memberships,candidate$memberships))/max(1,length(union(old$memberships,candidate$memberships)))>=.8)return()
    picked[[length(picked)+1L]] <<- candidate
  }
  for(catalog in unique(vapply(candidates,function(t)t$catalog,character(1))))add(candidates[[which(vapply(candidates,function(t)t$catalog==catalog,logical(1)))[1]]])
  for(candidate in candidates)add(candidate)
  for(term in picked) {
    nodes <- rbind(nodes,data.frame(id=term$id,kind='function',feature=term$feature,label=term$label,view='Function',column=NA_integer_,weight=NA_real_,factor=NA_integer_,factor_r2=NA_real_,importance=0,catalog=term$catalog,q_value=term$q_value,p_value=term$p_value,status=if(term$q_value<.05)'BH-supported' else 'Annotation match',memberships=paste(term$memberships,collapse=';'),views=paste(term$views,collapse=';'),source_url=term$source_url,overlap=length(term$memberships)))
    for(feature in term$memberships)edges <- rbind(edges,data.frame(source=feature,target=term$id,value=1,p_value=NA_real_,n=NA_integer_,kind='annotation',view_pair='',q_value=NA_real_,gene='',catalog=term$catalog))
  }
  nodes <- network_layout(nodes,edges,names(lab$prepared))
  list(nodes=nodes,edges=edges,correlation_tests=correlation_tests,functional_tests=functional_tests,functions=picked,
    options=list(factor=factor,top_per_view=top_per_view,threshold=threshold,q_max=q_max,term_limit=term_limit,include_within=include_within,include_identity=include_identity))
}
plot_integrated_network <- function(net,labels='key') {
  nodes <- net$nodes;edges <- net$edges;edges$x <- nodes$x[match(edges$source,nodes$id)];edges$y <- nodes$y[match(edges$source,nodes$id)];edges$xend <- nodes$x[match(edges$target,nodes$id)];edges$yend <- nodes$y[match(edges$target,nodes$id)]
  edges$evidence <- ifelse(edges$kind=='correlation',ifelse(edges$value>=0,'Positive correlation','Negative correlation'),ifelse(edges$kind=='annotation','Function membership','Same gene ID'))
  nodes$group <- ifelse(nodes$kind=='function',ifelse(nodes$q_value<.05,'BH-supported function','Annotation match'),nodes$view)
  nodes$size <- ifelse(nodes$kind=='function',4.5,2.2+1.4*sqrt(nodes$importance))
  counts <- ave(seq_len(nrow(nodes)),nodes$view,FUN=seq_along);show <- nodes$kind=='function'|labels=='all'|counts<=6
  ggplot2::ggplot()+ggplot2::geom_curve(data=edges,ggplot2::aes(x=x,y=y,xend=xend,yend=yend,color=evidence,linetype=evidence),curvature=.05,linewidth=.35,alpha=.38)+
    ggplot2::scale_color_manual(values=c('Positive correlation'='#42968e','Negative correlation'='#c07868','Function membership'='#a5a7b3','Same gene ID'='#8a9da5'))+
    ggplot2::scale_linetype_manual(values=c('Positive correlation'='solid','Negative correlation'='solid','Function membership'='dotted','Same gene ID'='dashed'))+
    ggplot2::geom_point(data=nodes,ggplot2::aes(x,y,fill=group,size=size,shape=kind),color='white',stroke=.7)+ggplot2::scale_shape_manual(values=c(feature=21,`function`=22))+
    ggplot2::scale_size_identity()+ggplot2::scale_fill_manual(values=c(view_colors,'BH-supported function'='#4f8c7c','Annotation match'='#a2a8bc'))+
    ggrepel::geom_text_repel(data=nodes[show,,drop=FALSE],ggplot2::aes(x,y,label=label),seed=557,size=2.8,max.overlaps=Inf,min.segment.length=0)+
    ggplot2::labs(title='Multi-omics functional network',subtitle=paste(sum(nodes$kind=='feature'),'features |',sum(nodes$kind=='function'),'function nodes |',sum(edges$kind=='correlation'),'correlation links'),color='Link evidence',fill='Node type',shape='Node')+
    ggplot2::theme_void()+ggplot2::theme(legend.position='bottom',legend.box='vertical',plot.title=ggplot2::element_text(face='bold',color='#111820',size=15),plot.subtitle=ggplot2::element_text(color='#607b83',margin=ggplot2::margin(b=12)),plot.margin=ggplot2::margin(20,20,20,20))
}
