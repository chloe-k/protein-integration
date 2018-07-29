# integrative DRW on combined feature data (updated in 2018/07/28)

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.


# edge direction (GMR_3_d) (202.30.3.222 server)
# m -> g (only anticorr)
# p -> g


# Classifier : rf(Random Forest)

################################## Result 35 in GMR_d ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = 10)

id_list <- c("35_1", "35_2", "35_3", "35_4",
             "35_5", "35_6", "35_7", "35_8",
             "35_9", "35_10", "35_11", "35_12",
             "35_13", "35_14", "35_15", "35_16")

prob_list <- c(0.2, 0.2, 0.2, 0.2,
               0.4, 0.4, 0.4, 0.4,
               0.6, 0.6, 0.6, 0.6,
               0.8, 0.8, 0.8, 0.8)

Gamma_list <- c(0.2, 0.4, 0.6, 0.8,
                0.2, 0.4, 0.6, 0.8,
                0.2, 0.4, 0.6, 0.8,
                0.2, 0.4, 0.6, 0.8)

make_GMR_d_model(id=id_list[1], prob = prob_list[1], Gamma = Gamma_list[1], mode = "GMR_3_d", AntiCorr = TRUE)
pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_35 <- foreach(i=2:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], prob = prob_list[i], Gamma = Gamma_list[i], mode = "GMR_3_d", AntiCorr = TRUE)
}

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_35_', i, '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}


for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR_d'), collapse = '')
  write.SigFeatures(res_fit=res_models[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}


# Plot for GMR model
title <- c("Result 35 GMR_d")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.02, perf_max = perf_max+0.02)
perf_heatmap(title, res_models, prob_list = prob_list, Gamma_list = Gamma_list)
