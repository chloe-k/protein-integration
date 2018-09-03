plotPerf <-
  function(xlabs, perf, title, measure, perf_min, perf_max, color="Dark2", group){
    library(scales)
    method = rownames(perf)
    pf <- as.matrix(perf)
    pf.means <- apply(pf,1,mean)
    pf.sd <- apply(pf,1,sd)
    pf.n <- apply(pf,1,length)
    pf.se <- pf.sd/sqrt(pf.n)

    print(list(pf.means, pf.sd, pf.se))
    
    data_AUC <- data.frame(Method=method,
                           meanscore=pf.means,
                           sd=pf.sd,
                           n=pf.n,
                           se=pf.se)
    
    data_AUC$Method_ <- factor(data_AUC$Method, as.character(data_AUC$Method))
    # data_AUC$group <- factor(data_AUC$group, levels = c("{p(G),p(M),p(P)}", "{p(GM)}",

    # model_name <- c(expression(iDRW(~G^R~M^R)),
    #                 expression(iDRW(~G^R~M^R~P)),
    #                 expression(iDRW_prop(~G^R~M^R~P)))
    
    # xlabs <- c("DRW(G+P)_mean", "DRW(G+P)_median", "DRW(G+P)_concat",
    #            "iDRW(GM)", "iDRW(GMP)", "iDRW-D^p(GMP)")
    
    # model_name <- c("DRW(G+P)_mean",
    #                 "DRW(G+P)_median",
    #                 "DRW(G+P)_concat",
    #                 "iDRW(GM)",
    #                 "iDRW(GMP)",
    #                 expression(iDRW(GMP)_~prop^(p)~))
    
    
    g <- ggplot(data=data_AUC, aes(x=Method_, y=meanscore, fill=factor(Method_, as.character(Method_)))) +
      geom_bar(aes(x=Method_), data=data_AUC, stat='identity', position=position_dodge()) +
      xlab('Model')+
      ylab(measure) +
      scale_color_brewer(palette=color) +
      scale_y_continuous(limits=c(perf_min,perf_max), oob=rescale_none)+
      theme(axis.text.x=element_text(angle=45, hjust=1, size=12), legend.position="none") +
      # geom_errorbar(aes(ymin=meanscore-se, ymax=meanscore+se), colour="black", width=.2, position=position_dodge(0.9)) +
      ggtitle(title)
    
    # g <- ggplot(data=data_AUC, aes(x=Method_, y=meanscore, fill=group)) +
    # g <- ggplot(data=df, aes(x=model, y=Accuracy, fill=factor(model, as.character(model)))) +
      
      # geom_bar(aes(x=model), data=Accuracy, stat='identity', position=position_dodge()) +
      
      # xlab('Model')+
      # ylab('Accuracy') +
      # scale_color_brewer(palette=color) +
      # # scale_fill_discrete("Set of\npathway profiles\nused in model") +
      # scale_y_continuous(limits=c(perf_min-0.01, perf_max+0.01), oob=rescale_none)+
      # theme(axis.text.x=element_text(angle=45, hjust=1, size=10), legend.title = element_text(hjust=0.5, size=11)) +
      # # scale_x_discrete(labels=model_name)+
      # # geom_errorbar(aes(ymin=meanscore-se, ymax=meanscore+se), colour="black", width=.2, position=position_dodge(0.9)) +
      # ggtitle(title)
    
    return(g)
  }

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}