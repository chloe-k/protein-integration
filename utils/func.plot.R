perf_boxplot <- function(title, xlabs, res_models, perf_min, perf_max, baseline=NULL) {
  df_list = list()
  for(i in 1:length(xlabs)) {
    #df_list[[i]] = data.frame(model=xlabs[i], Accuracy=res_models[[i]]$resample$Accuracy)
    df_list[[i]] = data.frame(model=xlabs[i], Accuracy=res_models[[i]]$results$Accuracy)
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
    # scale_x_continuous(name = 'restart prob(p-prob, g-Gamma)', 
    #                    breaks = c('p=0.2', 'p=0.2', 'p=0.2', 'p=0.2', 'p=0.4', 'p=0.4', 'p=0.4', 'p=0.4', 'p=0.6', 'p=0.6', 'p=0.6', 'p=0.6', 'p=0.8', 'p=0.8', 'p=0.8', 'p=0.8'), 
    #                    labels = c('p=0.2\ng=0.2', 'p=0.2\ng=0.4', 'p=0.2\ng=0.6', 'p=0.2\ng=0.8', 'p=0.4\ng=0.2', 'p=0.4\ng=0.4', 'p=0.4\ng=0.6', 'p=0.4\ng=0.8',
    #                               'p=0.6\ng=0.2', 'p=0.6\ng=0.4', 'p=0.6\ng=0.6', 'p=0.6\ng=0.8', 'p=0.8\ng=0.2', 'p=0.8\ng=0.4', 'p=0.8\ng=0.6', 'p=0.8\ng=0.8'))
  
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