################################## Result 28 in GR ############################################################
registerDoParallel(cores = 4)

# make RData after DRW

id_list <- c("28_0.2", "28_0.4", "28_0.6", "28_0.8")
Gamma_list <- c(0.2, 0.4, 0.6, 0.8, 0.9)


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gr_28 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GR_model(id=id_list[i], prob = 0.001, Gamma = Gamma_list[i])
}

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GR_', id_list[i], '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}

for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GR'), collapse = '')
  write.SigFeatures(res_fit=res_models[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}
