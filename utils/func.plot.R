perf_boxplot <- function(title, xlabs, res_models, perf_min, perf_max, baseline=NULL) {
  df_list = list()
  for(i in 1:length(xlabs)) {
    # df_list[[i]] = data.frame(model=xlabs[i], Accuracy=res_models[[i]]$resample$Accuracy)
    df_list[[i]] = data.frame(model=xlabs[i], Accuracy=max(res_models[[i]]$results$Accuracy))
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
    # geom_hline(aes_string(yintercept=baseline), linetype="dashed") +
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

perf_lineplot <- function(title, xlabs, res_models, perf_min, perf_max, Gamma_list) {
  
  df_list = list()
  for(i in 1:length(res_models)) {
    df_list[[i]] = data.frame(model=xlabs[i], Accuracy=mean(res_models[[i]]$resample$Accuracy), Gamma = Gamma_list[i], sd=sd(res_models[[i]]$resample$Accuracy))
    # df_list[[i]] = data.frame(model=xlabs[i], Accuracy=max(res_models[[i]]$results$Accuracy))
  }
  df_tot = Reduce(rbind, df_list)
  df <- df_tot[order(df_tot$Gamma, -df_tot$Accuracy),]
  model_order <- unique(df$model)
  df$model <- factor(df$model, levels = model_order)
  
  p <- ggplot(df, aes(x=Gamma, y=Accuracy, group=model, color=model)) +
    geom_line() +
    geom_point() +
    theme(plot.title = element_text(hjust = 0.5)) 
  print(p + ggtitle(title))
}

perf_lineplot_d <- function(title, xlabs, res_models, perf_min, perf_max, prob_list) {
  
  df_list = list()
  for(i in 1:length(res_models)) {
    df_list[[i]] = data.frame(model=xlabs[i], Accuracy=mean(res_models[[i]]$resample$Accuracy), Prob = prob_list[i], sd=sd(res_models[[i]]$resample$Accuracy))
    # df_list[[i]] = data.frame(model=xlabs[i], Accuracy=max(res_models[[i]]$results$Accuracy))
  }
  df_tot = Reduce(rbind, df_list)
  df <- df_tot[order(df_tot$Prob, -df_tot$Accuracy),]
  model_order <- unique(df$model)
  df$model <- factor(df$model, levels = model_order)
  
  p <- ggplot(df, aes(x=Prob, y=Accuracy, group=model, color=model)) +
    geom_line() +
    geom_point() +
    theme(plot.title = element_text(hjust = 0.5)) 
  print(p + ggtitle(title))
}

perf_facet_boxplot <- function(title, xlabs, res_models, perf_min, perf_max, baseline=NULL, prob_list, Gamma_list) {
  df_list = list()
  for(i in 1:length(xlabs)) {
    df_list[[i]] = data.frame(Accuracy=max(res_models[[i]]$results$Accuracy), P=prob_list[i], Gamma=Gamma_list[i])
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

perf_heatmap <- function(title, res_models, prob_list, Gamma_list){
  df_list = list()
  for(i in 1:length(res_models)) {
    df_list[[i]] = data.frame(Accuracy=max(res_models[[i]]$results$Accuracy), prob = prob_list[i], Gamma = Gamma_list[i])
    # df_list[[i]] = data.frame(model=xlabs[i], Accuracy=mean(res_models[[i]]$resample$Accuracy))
    # df_list[[i]] = data.frame(model=xlabs[i], Accuracy=max(res_models[[i]]$results$Accuracy))
  }
  df = Reduce(rbind, df_list)
  mat <- matrix(df$Accuracy, length(unique(df$prob)), length(unique(df$Gamma)), byrow = TRUE)
  rownames(mat) <- sprintf("p = %s", as.character(c(0.001, 0.01, 0.2, 0.4, 0.6, 0.8)))
  colnames(mat) <- sprintf("g = %s", as.character(c(0, 0.2, 0.4, 0.6, 0.8)))
  mat_breaks <- seq(min(mat), max(mat), length.out = 10) 
  
  pheatmap(mat, display_numbers = TRUE, number_format = "%.3f",
           main = title, cluster_cols = FALSE, cluster_rows = FALSE,
           legend = TRUE, cellwidth = 50, cellheight = 50, breaks = mat_breaks,
           color =  colorRampPalette(c("yellow", "red"))(10))
           # color = inferno(16))
  
  
}