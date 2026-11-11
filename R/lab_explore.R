# R equivalents of the browser's source-backed exploratory calculations.
factor_enrichment<-function(lab,mofa,annotations,factor=1,view='RNA',catalog='GO_BP',top_n=20,direction='absolute',mode='ora'){
 if((view=='Metabolomics')!=(catalog=='Metabolite_pathways'))stop('Genes use GO/Reactome; metabolites use compound pathways.')
 terms<-annotations$catalogs[[catalog]];if(is.null(terms))stop('Unknown catalog.')
 members<-lapply(terms,function(t)unlist(t$members,use.names=FALSE));annotated<-unique(unlist(members,use.names=FALSE))
 ids<-colnames(lab$prepared[[view]]$model);w<-mofa$weights[[view]][,factor];prefix<-paste0(view,'::');names(w)<-sub(prefix,'',names(w),fixed=TRUE)
 if(direction=='positive')w<-w[w>0] else if(direction=='negative')w<-w[w<0]
 ranked<-if(direction=='negative')sort(w) else if(direction=='positive')sort(w,decreasing=TRUE) else w[order(-abs(w))]
 selected<-names(head(ranked,top_n));universe<-intersect(ids,annotated);query<-intersect(selected,universe)
 if(mode=='coverage'){
  available<-intersect(colnames(lab$prepared[[view]]$retained),annotated)
  rows<-lapply(seq_along(terms),function(i){hit<-intersect(members[[i]],available);data.frame(id=terms[[i]]$id,name=terms[[i]]$name,overlap=length(hit),features=paste(sort(hit),collapse=';'))})
  d<-do.call(rbind,rows);d<-d[d$overlap>0,];d<-d[order(-d$overlap,d$name),];return(list(mode=mode,rows=d,measured_mapped=length(available),measured_total=ncol(lab$prepared[[view]]$retained)))
 }
 N<-length(universe);n<-length(query);min_size<-if(view=='Metabolomics')2 else 3
 rows<-lapply(seq_along(terms),function(i){set<-intersect(members[[i]],universe);M<-length(set);if(M<min_size)return(NULL);hit<-intersect(set,query);k<-length(hit);data.frame(id=terms[[i]]$id,name=terms[[i]]$name,overlap=k,query_size=n,background_size=N,term_in_background=M,p_value=if(n)phyper(k-1,M,N-M,n,lower.tail=FALSE) else 1,features=paste(sort(hit),collapse=';'))})
 d<-do.call(rbind,rows);if(is.null(d))d<-data.frame(id=character(),name=character(),overlap=integer(),query_size=integer(),background_size=integer(),term_in_background=integer(),p_value=double(),features=character())
 d$q_value<-p.adjust(d$p_value,'BH');d$gene_ratio<-if(n)d$overlap/n else 0;d<-d[order(d$q_value,d$p_value,-d$overlap),]
 list(mode='ora',rows=d,universe_size=N,query_mapped=n,query_selected=length(selected),tested_terms=nrow(d),significant_terms=sum(d$q_value<.05),view=view,catalog=catalog,factor=factor)
}
plot_pathway_results<-function(e){
 d<-head(e$rows[e$rows$overlap>0,,drop=FALSE],15)
 if(!nrow(d))return(ggplot2::ggplot()+ggplot2::annotate('text',x=0,y=0,label='No mapped overlaps for this selection')+ggplot2::theme_void())
 d$name<-factor(d$name,levels=rev(d$name))
 if(e$mode=='coverage')p<-ggplot2::ggplot(d,ggplot2::aes(overlap,name,size=overlap))+ggplot2::geom_point(color='#C38B3B',alpha=.85)+ggplot2::labs(title='Mapped measured metabolites by pathway',subtitle='Coverage, not enrichment or pathway activity',x='Mapped measured compounds',y=NULL,size='Count')
 else p<-ggplot2::ggplot(d,ggplot2::aes(gene_ratio,name,size=overlap,color=-log10(pmax(q_value,1e-300))))+ggplot2::geom_point(alpha=.85)+ggplot2::scale_color_gradient(low='#92BAB0',high='#566E9E')+ggplot2::labs(title='Top-weight gene-set or pathway candidates',subtitle=paste(e$view,'| BH across',e$tested_terms,'tested terms |',e$significant_terms,'q < 0.05'),x='Fraction of mapped query in term',y=NULL,size='Overlap',color='-log10(q)')
 p+ggplot2::scale_size_area(max_size=8)+lab_theme()+ggplot2::theme(axis.text.y=ggplot2::element_text(size=9),legend.position='right')
}
factor_network<-function(lab,mofa,factor=1,top_per_view=6,mode='factor',threshold=.6){
 nodes<-do.call(rbind,lapply(seq_along(lab$prepared),function(i){view<-names(lab$prepared)[i];w<-mofa$weights[[view]][,factor];ix<-head(order(-abs(w)),top_per_view);theta<-(-90+(i-1)*360/length(lab$prepared)-42+84*(seq_along(ix)-.5)/length(ix))*pi/180;data.frame(id=paste(view,colnames(lab$prepared[[view]]$model)[ix],sep='::'),feature=colnames(lab$prepared[[view]]$model)[ix],view=view,column=ix,weight=w[ix],x=cos(theta),y=sin(theta))}))
 edges<-data.frame(source=character(),target=character(),value=double())
 if(mode=='factor')edges<-data.frame(source=paste('Factor',factor),target=nodes$id,value=nodes$weight)
 else for(i in seq_len(nrow(nodes)-1))for(j in (i+1):nrow(nodes)){if(nodes$view[i]==nodes$view[j])next;x<-mofa$model_input[[nodes$view[i]]][,nodes$column[i]];y<-mofa$model_input[[nodes$view[j]]][,nodes$column[j]];ok<-is.finite(x)&is.finite(y);if(sum(ok)<10)next;r<-cor(x[ok],y[ok],method='spearman');if(is.finite(r)&&abs(r)>=threshold)edges<-rbind(edges,data.frame(source=nodes$id[i],target=nodes$id[j],value=r))}
 if(mode=='correlation'&&nrow(edges)>100)edges<-edges[head(order(-abs(edges$value)),100),]
 list(nodes=nodes,edges=edges,factor=factor,mode=mode)
}
plot_factor_network<-function(net){
 nodes<-net$nodes;pos<-rbind(nodes[,c('id','x','y')],data.frame(id=paste('Factor',net$factor),x=0,y=0));edges<-net$edges
 edges$x<-pos$x[match(edges$source,pos$id)];edges$y<-pos$y[match(edges$source,pos$id)];edges$xend<-pos$x[match(edges$target,pos$id)];edges$yend<-pos$y[match(edges$target,pos$id)];edges$sign<-ifelse(edges$value>=0,'Positive','Negative')
 p<-ggplot2::ggplot()+ggplot2::geom_curve(data=edges,ggplot2::aes(x=x,y=y,xend=xend,yend=yend,color=sign),curvature=.12,alpha=.4,linewidth=.6)+ggplot2::scale_color_manual(values=c(Positive='#579A99',Negative='#C78265'))+ggplot2::geom_point(data=nodes,ggplot2::aes(x,y,fill=view),shape=21,size=5,color='white',stroke=.8)+ggplot2::scale_fill_manual(values=view_colors)+ggrepel::geom_text_repel(data=nodes,ggplot2::aes(x,y,label=feature),max.overlaps=Inf,size=3,min.segment.length=0,seed=557)+ggplot2::coord_equal(xlim=c(-1.65,1.65),ylim=c(-1.3,1.3))+ggplot2::labs(title=paste('Factor',net$factor,'across molecular layers'),subtitle=if(net$mode=='factor')'Edges are model weights; feature selection is within each view' else 'Cross-view observed Spearman correlations; not physical interactions')+lab_theme()+ggplot2::theme(axis.title=ggplot2::element_blank(),axis.text=ggplot2::element_blank(),panel.grid=ggplot2::element_blank(),panel.grid.major=ggplot2::element_blank(),panel.grid.minor=ggplot2::element_blank(),legend.position='bottom')
 if(net$mode=='factor')p<-p+ggplot2::annotate('point',x=0,y=0,size=17,color='#173F4C')+ggplot2::annotate('text',x=0,y=0,label=paste('F',net$factor,sep=''),color='white',size=4)
 p
}
compare_factor_counts<-function(lab,counts=c(2,3,5),scale_views=FALSE,runtime='managed'){
 rows<-lapply(counts,function(k){m<-fit_lab_mofa(lab,factors=k,scale_views=scale_views,runtime=runtime,outfile=tempfile(fileext='.hdf5'));data.frame(K=k,view=names(m$total_r2),total_r2=as.numeric(m$total_r2))});do.call(rbind,rows)
}
plot_factor_count_comparison<-function(d){ggplot2::ggplot(d,ggplot2::aes(K,total_r2,color=view))+ggplot2::geom_line(linewidth=.8)+ggplot2::geom_point(size=3)+ggplot2::scale_color_manual(values=view_colors)+ggplot2::scale_y_continuous(labels=scales::percent)+ggplot2::labs(title='Reconstruction across actual factor-count refits',subtitle='Same features and preparation; in-sample fit, not an optimal-K rule',x='Number of factors K',y='Total reconstruction R-squared')+lab_theme()}
