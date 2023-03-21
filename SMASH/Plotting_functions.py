pattern_plot_SMASH <- function(pltdat, gene_name, xy=T,main=F,titlesize=2,pointsize=3,min.pand=0.99,max.pand=1.01,title=NULL,pal=NULL,expand_par=0.05,ncolors=5){
  if(!xy){
    xy              <- matrix(as.numeric(do.call(rbind,strsplit(as.character(pltdat[,1]),split="x"))),ncol=2)
    rownames(xy)    <- as.character(pltdat[,1])
    colnames(xy)    <- c("x","y")
    pd              <- cbind.data.frame(xy,pltdat[,2:ncol(pltdat)])
  }else{
    pd              <- pltdat
  }
  library(viridis)
  if(is.null(pal)){
    pal             <- colorRampPalette(c("antiquewhite",viridis_pal()(10)[c(10, 9, 8, 7)]))
    }
  
    gpt             <- ggplot(pd,aes(x=x, y=y,color=pd[,3])) + 
    geom_point(size=pointsize) +
    scale_color_gradientn(colours=pal(ncolors), name = "Relative expression")+ 
    scale_x_discrete(expand = c(0, expand_par))+ 
    scale_y_discrete(expand = c(0, expand_par))+
    expand_limits(x=c(min(pd$x)*min.pand,max(pd$x)*max.pand),y=c(min(pd$y)*min.pand,max(pd$y)*max.pand))+ 
    theme_bw()
  if(main){
    if(is.null(title)){
      title=colnames(pd)[igene+2]
    }
    out = gpt + labs(title = title, x = NULL, y = NULL)+ scale_y_reverse() + 
      theme(legend.position = "bottom", plot.title = element_text(hjust = 0.5,size=rel(titlesize),face="italic"), 
            axis.line=element_blank(),axis.text.x=element_blank(), axis.text.y=element_blank(),
            axis.ticks.x=element_blank(),axis.ticks.y=element_blank(), 
            plot.title = element_text(hjust = 0.5))+ ggtitle(gene_name)
  }else{
    out = gpt + labs(title = NULL, x = NULL, y = NULL) + scale_y_reverse() +  scale_fill_discrete(name = "Value")+
      theme(legend.position = "bottom", axis.line=element_blank(),axis.text.x=element_blank(), 
            axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank(), 
            plot.title = element_text(hjust = 0.5)) + ggtitle(gene_name)
  }
  return(out)
}