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

