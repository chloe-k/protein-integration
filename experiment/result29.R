################################## Result 29 in GR_d ############################################################
registerDoParallel(cores = 4)
id_list <- c("29_1", "29_2", "29_3", "29_4",
             "29_5", "29_6", "29_7", "29_8",
             "29_9", "29_10", "29_11", "29_12",
             "29_13", "29_14", "29_15", "29_16")


prob_list <- c(0.2, 0.2, 0.2, 0.2,
               0.4, 0.4, 0.4, 0.4,
               0.6, 0.6, 0.6, 0.6,
               0.8, 0.8, 0.8, 0.8)

Gamma_list <- c(0.2, 0.4, 0.6, 0.8,
                0.2, 0.4, 0.6, 0.8,
                0.2, 0.4, 0.6, 0.8,
                0.2, 0.4, 0.6, 0.8)


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gr_d <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GR_d_model(id=id_list[i], prob = prob_list[i], Gamma = Gamma_list[i])
}

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GR_d_', id_list[i], '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}

result_name <- paste(c('result',id_list[1],'_GR_d'), collapse = '')
write.SigFeatures(res_fit=res_models[[1]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)

result_name <- paste(c('result',id_list[16],'_GR_d'), collapse = '')
write.SigFeatures(res_fit=res_models[[16]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)


for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GR_d'), collapse = '')
  write.SigFeatures(res_fit=res_models[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

title <- c("Result 29 GR_d")

perf_heatmap(title, res_models, prob_list = prob_list, Gamma_list = Gamma_list)
