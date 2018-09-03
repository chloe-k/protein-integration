# integrative DRW on combined feature data (updated in 2018/07/28)

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.


# edge direction (GMR_3)
# m -> g (only anticorr)
# p -> g


# Classifier : rf(Random Forest)


################################## Result 34 in GMR ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = 10)

id_list <- c("34_1", "34_2", "34_3", "34_4")
Gamma_list <- c(0.2, 0.4, 0.6, 0.8)

make_GMR_model(id=id_list[1], prob = 0.001, Gamma = Gamma_list[1], mode = "GMR_3", AntiCorr = TRUE)
pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_34 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], prob = 0.001, Gamma = Gamma_list[i], mode = "GMR_3", AntiCorr = TRUE)
}

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_34_', i, '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}


for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR'), collapse = '')
  write.SigFeatures(res_fit=res_gmr[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}


# Plot for GMR model
title <- c("Result 34_GMR")
xlabs <- c("g=0.2", "g=0.4", "g=0.6", "g=0.8")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.02, perf_max = perf_max+0.02)
perf_lineplot(title, xlabs, res_models, perf_min, perf_max, Gamma_list)



################################## Result 34 in GMR_d ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = 10)

id_list <- c("34_1_d", "34_2_d", "34_3_d", "34_4_d",
             "34_5_d", "34_6_d", "34_7_d", "34_8_d",
             "34_9_d", "34_10_d", "34_11_d", "34_12_d",
             "34_13_d", "34_14_d", "34_15_d", "34_16_d")

prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), each = 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)

make_GMR_d_model(id=id_list[1], prob = prob_list[1], Gamma = Gamma_list[1], mode = "GMR_3_d", AntiCorr = TRUE)
pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_34_d <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], prob = prob_list[i], Gamma = Gamma_list[i], mode = "GMR_3_d", AntiCorr = TRUE)
}

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}


for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR_d'), collapse = '')
  write.SigFeatures(res_fit=res_gmr[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}


# Plot for GMR model
title <- c("Result 34_GMR_d")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_heatmap(title, res_models, prob_list = prob_list, Gamma_list = Gamma_list)

