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


res_gmr_33 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
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
