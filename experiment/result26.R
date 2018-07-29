# integrative DRW on combined feature data (updated in 2018/07/20)
# concat directed pathway graphs within each profile (GMR)

# All gene symbols are converted to Entrez gene id
# 5-fold CV(10 iters) was performed for tuning parameter in Random Forest.
# 5-fold CV(10 iters) was performed for get top N pathways.
# LOOCV was performed for model evaluation

# Dppigraph(Entrez).rda was used 
# Dup_rppa_data.RData was used in GM, GMR -> unmapped genes with rppa proteins are removed in rnaseq, imputed_methyl profiles
# Entrez_data.RData was used in GM_base


# data dimension
# rnaseq - 183 genes, 376 samples
# imputed_methyl - 176 genes, 376 samples
# rppa - 184 genes, 376 samples

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)




########################################### Result26 #######################################
# res_models <- list()
# #------------------------- RNAseq + Methyl -------------------------#
# gm <- g %du% m
# testStatistic <- c("DESeq2", "t-test")
# profile_name <- c("rna(res26)", "meth(res26)")
# x=list(rnaseq, imputed_methyl)
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result26_GM", prob = 0.001, Gamma = 0.4, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)
# 
# res_pa_GM_26 <- fit.classification(y=y, samples = samples, id = "result26_GM", datapath = datapath, respath = respath,
#                                    profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                                    nFolds = 5, numTops=50, iter = 10)
# 
# save(res_pa_GM_26, file=file.path('data/model/res_pa_GM_26.RData'))
# 
# write.SigFeatures(res_fit=res_pa_GM_26, id="result26_GM", profile_name=profile_name, method="DRW", respath=respath)
# 
# res_models <- c(res_models, list(res_pa_GM_26))
# 
# #------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
# gmr <- g %du% m %du% r
# testStatistic <- c("DESeq2", "t-test", "t-test")
# profile_name <- c("rna(res26)", "meth(res26)", "rppa(res26)")
# x=list(rnaseq, imputed_methyl, rppa)
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result26_GMR", prob = 0.001, Gamma = 0.9, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)
# 
# res_pa_GMR_26 <- fit.classification(y=y, samples = samples, id = "result26_GMR", datapath = datapath, respath = respath,
#                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                                     nFolds = 5, numTops=50, iter = 10)
# 
# save(res_pa_GMR_26, file=file.path('data/model/res_pa_GMR_26.RData'))
# 
# write.SigFeatures(res_fit=res_pa_GMR_26, id="result26_GMR", profile_name=profile_name, method="DRW", respath=respath)
# 
# res_models <- c(res_models, list(res_pa_GMR_26))
# 
# #------------------------- RNAseq + Methyl (baseline) -------------------------#
# gm <- g %du% m
# testStatistic <- c("DESeq2", "t-test")
# profile_name <- c("rna(Entrez)", "meth(Entrez)")
# x=list(rnaseq, imputed_methyl)
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result26_base_GM", prob = 0.001, Gamma = 0.4, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)
# 
# res_pa_GM_26_base <- fit.classification(y=y, samples = samples, id = "result26_base_GM", datapath = datapath, respath = respath,
#                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                                         nFolds = 5, numTops=50, iter = 10)
# 
# save(res_pa_GM_26_base, file=file.path('data/model/res_pa_GM_26_base.RData'))
# 
# write.SigFeatures(res_fit=res_pa_GM_26_base, id="result26_base_GM", profile_name=profile_name, method="DRW", respath=respath)
# 
# baseline <- max(res_pa_GM_26_base$results$Accuracy)
# 
# # Plot for Result26
# 
# title <- c("Result 26")
# xlabs <- c("GM", "GMR")
# 
# perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
# perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.05, perf_max = perf_max+0.05, baseline)
# 


################################## Result 26 in GMR ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = 10)

id_list <- c("26_0.2", "26_0.4", "26_0.6", "26_0.8", "26_0.9")
Gamma_list <- c(0.2, 0.4, 0.6, 0.8, 0.9)

make_GMR_model(id=id_list[1], prob = 0.001, Gamma = Gamma_list[1], mode = "GMR_26")

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

foreach(i=2:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], prob = 0.001, Gamma = Gamma_list[i], mode = "GMR_26")
}

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_', id_list[i], '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}


title <- c("Result 26 GMR")
xlabs <- c("[g=0.2]", "[g=0.4]", "[g=0.6]", "[g=0.8]", "[g=0.9]")
perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_lineplot(title = title, xlabs = xlabs, res_models = res_models, perf_max = perf_max, perf_min = perf_min, Gamma_list = Gamma_list)



################################## Result 26 in GMR_d ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = 6)

id_list <- c("26_1", "26_2", "26_3", "26_4",
             "26_5", "26_6", "26_7", "26_8", 
             "26_9", "26_10", "26_11", "26_12",
             "26_13", "26_14", "26_15", "26_16")

prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), each = 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)

make_GMR_d_model(id=id_list[1], prob = prob_list[1], Gamma = Gamma_list[1], mode = "GMR_d_26")

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_d_26 <- foreach(i=2:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], prob = prob_list[i], Gamma = Gamma_list[i], mode = "GMR_d_26")
}


res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}


# Plot for GMR_d_26 
title <- c("Result 26 GMR_d Heatmap")
perf_heatmap(title, res_models, prob_list = prob_list, Gamma_list = Gamma_list)

# Plot for comparison with GMR_26, GMR_d_26, GM_baseline
title <- c("Result 26")
gmr <- get(load('data/model/res_pa_GMR_26_0.4_LOOCV.RData'))
gmr_d <- get(load('data/model/res_pa_GMR_d_26_13_LOOCV.RData'))
res_models <- list(gmr, gmr_d)

gm <- get(load('data/model/res_pa_GM_26_base.RData'))
baseline <- max(gm$results$Accuracy)

xlabs <- c("GMR", "GMR_d")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.01, perf_max = perf_max+0.01, baseline)
