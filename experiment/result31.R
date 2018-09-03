# integrative DRW on combined feature data (updated in 2018/07/28)
# concat directed pathway graphs within each profile (G & M & R & GM & GR & MP & GMR)

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.


# edge direction (GMR_1)
# p -> m
# p -> g

# Classifier : rf(Random Forest)


################################## Result 31 in GMR ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = num_cores)

id_list <- c("31_1", "31_2", "31_3", "31_4", "31_5", "31_6", "31_7")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_31 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.6)
}

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_31_', i, '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}

for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR'), collapse = '')
  write.SigFeatures(res_fit=res_gmr[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

# Plot for GMR model
title <- c("Result 31_GMR")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.02, perf_max = perf_max+0.02)

