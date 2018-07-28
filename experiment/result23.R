# integrative DRW on combined feature data (updated in 2018/07/24)
# concat directed pathway graphs within each profile (G & M & R & GM & GR & MP & GMR)

# g=0.9 was used for restart probability
# When pathway activity score was calculated, each type of weight was used.
# This experiment is based on GM or GMR_d model
# G -> perform GMR_d model and only 'g' type weight is included for pathwway activity score calculation 

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.


# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)



################################## Result 23 in GM ############################################################


################################## Result 23 in GMR ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = num_cores)

id_list <- c("23_1", "23_2", "23_3", "23_4", "23_5", "23_6", "23_7")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_23 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.9)
}


for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_23_', i, '_LOOCV.RData'), collapse = ''))
}

res_gmr <- list(res_pa_GMR_23_1_LOOCV, res_pa_GMR_23_2_LOOCV, res_pa_GMR_23_3_LOOCV, res_pa_GMR_23_4_LOOCV, 
                res_pa_GMR_23_5_LOOCV, res_pa_GMR_23_6_LOOCV, res_pa_GMR_23_7_LOOCV)

for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR'), collapse = '')
  write.SigFeatures(res_fit=res_gmr[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

# Plot for GMR model
title <- c("Result 23_GMR")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr, perf_min = perf_min-0.02, perf_max = perf_max+0.02)




################################## Result 23 in GMR_d ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = 4)

id_list <- c("23_G", "23_M", "23_P", "23_GM", "23_GP", "23_MP", "23_GMP")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_23 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.4, Gamma = 0.2)
}


for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '_LOOCV.RData'), collapse = ''))
}

res_gmr_d <- list(res_pa_GMR_d_23_G_LOOCV, res_pa_GMR_d_23_M_LOOCV, res_pa_GMR_d_23_P_LOOCV,
                  res_pa_GMR_d_23_GM_LOOCV, res_pa_GMR_d_23_GP_LOOCV, res_pa_GMR_d_23_MP_LOOCV,
                  res_pa_GMR_d_23_GMP_LOOCV)


for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR_d'), collapse = '')
  write.SigFeatures(res_fit=res_gmr_d[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

# Plot for GMR_d model
title <- c("Result 23_GMR_d")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr_d, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_d, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr_d, perf_min = perf_min-0.02, perf_max = perf_max+0.02)


#########################################################################################################################################
# plot

# GM
title <- c("Result 23_GM")
xlabs <- c("G", "M")

perf_min <- min(sapply(X = res_gm_25, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gm_25, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gm_25, perf_min = perf_min-0.02, perf_max = perf_max+0.02)

#GMR_d
title <- c("Result 23_GMR_d")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr_d_25, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_d_25, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr_d_25, perf_min = perf_min-0.02, perf_max = perf_max+0.02)
