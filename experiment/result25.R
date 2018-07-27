################################## Result 25 ############################################################
num_cores <- detectCores()/2
registerDoParallel(cores = 4)


id_list <- c("25_5", "25_10", "25_15")
lim_list <- c(5, 10, 15)


# pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
#           "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")
# res_gmr_25 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
#   make_GMR_model(id=id_list[i], lim=lim_list[i], type_used = "gmp", prob = 0.001, Gamma = 0.9)
# }
# 
# res_gm_25 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
#   make_GM_model(id=id_list[i], lim = lim_list[i], type_used = "gmp", prob = 0.001, Gamma = 0.4)
# }


for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_', id_list[i], '_LOOCV.RData'), collapse = ''))
}

for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GM_', id_list[i], '_LOOCV.RData'), collapse = ''))
}

res_gm_25 <- list(res_pa_GM_25_5_LOOCV, res_pa_GM_25_10_LOOCV, res_pa_GM_25_15_LOOCV)
res_gmr_25 <- list(res_pa_GMR_25_5_LOOCV, res_pa_GMR_25_10_LOOCV, res_pa_GMR_25_15_LOOCV)

for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GM'), collapse = '')
  write.SigFeatures(res_fit=res_gm_25[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR'), collapse = '')
  write.SigFeatures(res_fit=res_gmr_25[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}



#########################################################################################################################################
# Plot for Result 25

res_model_25_5 <- list(res_pa_GM_25_5_LOOCV, res_pa_GMR_25_5_LOOCV)
res_model_25_10 <- list(res_pa_GM_25_10_LOOCV, res_pa_GMR_25_10_LOOCV)
res_model_25_15 <- list(res_pa_GM_25_15_LOOCV, res_pa_GMR_25_15_LOOCV)

title <- c("Result 25_5")
xlabs <- c("GM", "GMR")
perf_min <- min(sapply(X = res_model_25_5, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_5, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_5, perf_min = perf_min-0.05, perf_max = perf_max+0.05)


title <- c("Result 25_10")
xlabs <- c("GM", "GMR")
perf_min <- min(sapply(X = res_model_25_10, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_10, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_10, perf_min = perf_min-0.05, perf_max = perf_max+0.05)


title <- c("Result 25_15")
xlabs <- c("GM", "GMR")
perf_min <- min(sapply(X = res_model_25_15, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_15, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_15, perf_min = perf_min-0.05, perf_max = perf_max+0.05)

