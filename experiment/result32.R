# integrative DRW on combined feature data (updated in 2018/07/28)
# concat directed pathway graphs within each profile (G & M & R & GM & GR & MP & GMR)

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.


# edge direction
# m <-> g
# p <-> g

# Classifier : rf(Random Forest)


################################## Result 32 in GMR ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = num_cores)

id_list <- c("32_1", "32_2", "32_3", "32_4", "32_5", "32_6", "32_7")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_32 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.2, Gamma = 0.4, mode = "GMR_bidir")
}


for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_32_', i, '_LOOCV.RData'), collapse = ''))
}

res_gmr <- list(res_pa_GMR_32_1_LOOCV, res_pa_GMR_32_2_LOOCV, res_pa_GMR_32_3_LOOCV, res_pa_GMR_32_4_LOOCV, 
                res_pa_GMR_32_5_LOOCV, res_pa_GMR_32_6_LOOCV, res_pa_GMR_32_7_LOOCV)

for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR'), collapse = '')
  write.SigFeatures(res_fit=res_gmr[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

# Plot for GMR model
title <- c("Result 32_GMR")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr, perf_min = perf_min-0.02, perf_max = perf_max+0.02)

