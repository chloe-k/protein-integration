################################## Result 25 ############################################################

num_cores <- detectCores()/3
num_cores <- floor(num_cores*2)
cl <- makeCluster(num_cores)
registerDoParallel(cl = cl)


id_list <- c("25_5", "25_10", "25_15")
lim_list <- c(5, 10, 15)


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")
res_gmr_25 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], lim=lim_list[i])
}

res_gmr_d_25 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], lim = lim_list[i])
}

stopCluster(cl)



#########################################################################################################################################
# Plot for Result 25

res_models <- list(res_pa_GMR_25_5_LOOCV, res_pa_GMR_d_25_5_LOOCV, 
                   res_pa_GMR_25_10_LOOCV, res_pa_GMR_d_25_10_LOOCV,
                   res_pa_GMR_25_15_LOOCV, res_pa_GMR_d_25_15_LOOCV)

title <- c("Result 25")
xlabs <- c("GMR_5", "GMR_d_5", "GMR_10", "GMR_d_10", "GMR_15", "GMR_d_15")
perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.15, perf_max = perf_max+0.15)