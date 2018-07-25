################################## Result 27 in GMR ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = 4)


type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

id_list <- c("27_1_0.2", "27_2_0.2", "27_3_0.2", "27_4_0.2", "27_5_0.2", "27_6_0.2", "27_7_0.2")
res_gmr_27_0.2 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.2)
}

id_list <- c("27_1_0.4", "27_2_0.4", "27_3_0.4", "27_4_0.4", "27_5_0.4", "27_6_0.4", "27_7_0.4")
res_gmr_27_0.4 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.4)
}

id_list <- c("27_1_0.6", "27_2_0.6", "27_3_0.6", "27_4_0.6", "27_5_0.6", "27_6_0.6", "27_7_0.6")
res_gmr_27_0.6 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.6)
}

id_list <- c("27_1_0.8", "27_2_0.8", "27_3_0.8", "27_4_0.8", "27_5_0.8", "27_6_0.8", "27_7_0.8")
res_gmr_27_0.8 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.8)
}

id_list <- c("27_1_0.9", "27_2_0.9", "27_3_0.9", "27_4_0.9", "27_5_0.9", "27_6_0.9", "27_7_0.9")
res_gmr_27_0.9 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.9)
}


for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_27_', i, '_LOOCV.RData'), collapse = ''))
}



res_gmr <- list(res_pa_GMR_27_1_LOOCV, res_pa_GMR_27_2_LOOCV, res_pa_GMR_27_3_LOOCV, res_pa_GMR_27_4_LOOCV, 
                res_pa_GMR_27_5_LOOCV, res_pa_GMR_27_6_LOOCV, res_pa_GMR_27_7_LOOCV)

for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR'), collapse = '')
  write.SigFeatures(res_fit=res_gmr[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

# Plot for GMR model
title <- c("Result 27_GMR")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr, perf_min = perf_min-0.02, perf_max = perf_max+0.02)




################################## Result 27 in GMR_d ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = 4)

id_list <- c("27_1", "27_2", "27_3", "27_4", "27_5", "27_6", "27_7")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.4, Gamma = 0.2)
}


for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_d_27_', i, '_LOOCV.RData'), collapse = ''))
}

res_gmr_d <- list(res_pa_GMR_d_27_1_LOOCV, res_pa_GMR_d_27_2_LOOCV, res_pa_GMR_d_27_3_LOOCV, res_pa_GMR_d_27_4_LOOCV, 
                  res_pa_GMR_d_27_5_LOOCV, res_pa_GMR_d_27_6_LOOCV, res_pa_GMR_d_27_7_LOOCV)


for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR'), collapse = '')
  write.SigFeatures(res_fit=res_gmr_d[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

# Plot for GMR_d model
title <- c("Result 27_GMR_d")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr_d, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_d, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr_d, perf_min = perf_min-0.02, perf_max = perf_max+0.02)