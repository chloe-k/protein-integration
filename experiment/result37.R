# integrative DRW on combined feature data (updated in 2018/07/30)

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.

# Heat diffusion was performed for pathway graph (RPPA) diffusion


# edge direction (GMR_d)
# m -> g 
# p -> g 


# Classifier : rf(Random Forest)


################################## Result 37 in GMR_d ############################################################

registerDoParallel(cores = 5)

id_list <- c("37_1", "37_2", "37_3", "37_4",
             "37_5", "37_6", "37_7", "37_8",
             "37_9", "37_10", "37_11", "37_12",
             "37_13", "37_14", "37_15", "37_16")

prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), each = 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)

make_GMR_d_model(id=id_list[1], prob = prob_list[1], Gamma = Gamma_list[1], mode = "GMR_d_heat", AntiCorr = FALSE)
pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_37 <- foreach(i=2:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], prob = prob_list[i], Gamma = Gamma_list[i], mode = "GMR_d_heat", AntiCorr = FALSE)
}

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_37_', i, '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}


for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR_d'), collapse = '')
  write.SigFeatures(res_fit=res_models[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}


# Plot for GMR_d model
title <- c("Result 37 GMR_d")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_heatmap(title, res_models, prob_list = prob_list, Gamma_list = Gamma_list)
