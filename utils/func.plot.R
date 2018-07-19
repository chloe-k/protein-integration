perf_boxplot <- function(title, xlabs, res_models, perf_min, perf_max, baseline=NULL) {
  df_list = list()
  for(i in 1:length(xlabs)) {
    df_list[[i]] = data.frame(model=xlabs[i], Accuracy=res_models[[i]]$resample$Accuracy)
    # df_list[[i]] = data.frame(model=xlabs[i], Accuracy=res_models[[i]]$results$Accuracy)
  }
  df = Reduce(rbind, df_list)
  
  p <- ggplot(df, aes(x=model, y=Accuracy, fill=model)) +
    geom_boxplot() +
    xlab('Model') +
    ylab('Accuracy') +
    scale_color_brewer(palette="Dark2") +
    scale_y_continuous(limits=c(perf_min,perf_max)) +
    # geom_jitter(alpha=0.4, size=0.6, position=position_jitter(width=0.1,height=0)) +
    theme(axis.text.x=element_text(angle=45, hjust=1, size=12), legend.position="none") +
    #geom_hline(aes_string(yintercept=baseline), linetype="dashed") +
    theme(plot.title = element_text(hjust = 0.5)) 
  print(p + ggtitle(title))
}

perf_barplot <- function(xlabs, res_models, perf_min, perf_max, baseline=NULL) {
  
  res <- Reduce(rbind,lapply(X = res_models, FUN=function(x) x$results$Accuracy))
  row.names(res) <- xlabs
  
  g1 <- plotPerf(res, title = "", measure='Accuracy', perf_min=perf_min, perf_max=perf_max, color="Dark2")
  if(!is.null(baseline)) g1 <- g1 + geom_hline(aes_string(yintercept=baseline), linetype="dashed")
  
  print(g1)
  
}

perf_lineplot <- function(fname_res, perf_min, perf_max) {
  
  res <- read.table(file = fname_res,header = T)
  
  p <- ggplot(data=res, aes(x=k, y=Accuracy, group=model, colour=model)) +
    geom_line() +
    scale_y_continuous(limits=c(perf_min,perf_max)) +
    scale_color_brewer(palette="Dark2") +
    scale_x_continuous(breaks=seq(5,100,by=10)) +
    geom_point(aes(shape=model), size=1.5)
  
  print(p)
  
}

perf_facet_boxplot <- function(title, xlabs, res_models, perf_min, perf_max, baseline=NULL) {
  df_list = list()
  for(i in 1:length(xlabs)) {
    p <- -1
    g <- -1
    if(i%%5 == 1) g <- "0"
    else if(i%%5 == 2) g <- "0.2"
    else if(i%%5 == 3) g <- "0.4"
    else if(i%%5 == 4) g <- "0.6"
    else if(i%%5 == 0) g <- "0.8"
    
    if((i-1)%/%5 == 0) p <- "0.001"
    else if((i-1)%/%5 == 1) p <- "0.01"
    else if((i-1)%/%5 == 2) p <- "0.2"
    else if((i-1)%/%5 == 3) p <- "0.4"
    else if((i-1)%/%5 == 4) p <- "0.6"
    else if((i-1)%/%5 == 5) p <- "0.8"
    
    df_list[[i]] = data.frame(Accuracy=res_models[[i]]$resample$Accuracy, P=p, Gamma=g)
  }
  df = Reduce(rbind, df_list)
  
  facet_p <- ggplot(df, aes(x=P, y=Accuracy, group=P)) +
    geom_boxplot() +
    xlab('Restart Prob in network diffusion(p)') +
    ylab('Accuracy') +
    scale_color_brewer(palette="Dark2") +
    scale_y_continuous(limits=c(perf_min,perf_max)) +
    #theme(axis.text.x=element_text(angle=45, hjust=1, size=12), legend.position="none") +
    geom_hline(aes_string(yintercept=baseline), linetype="dashed") +
    theme(plot.title = element_text(hjust = 0.5))
  facet_p <- facet_p + facet_grid(Gamma ~ .)
  print(facet_p + ggtitle(title))
}