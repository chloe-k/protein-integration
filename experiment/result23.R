# integrative DRW on combined feature data (updated in 2018/07/20)
# concat directed pathway graphs within each profile (G & M & R & GM & GR & MP & GMR)

# p=0.4, g=0.4 was used for restart probability
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


################################## Result 23 in GMR_d ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = num_cores)

id_list <- c("23_G", "23_M", "23_P", "23_GM", "23_GP", "23_MP", "23_GMP")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_23 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.4, Gamma = 0.2)
}



# #------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
# gmr <- list(g, m, r)
# testStatistic <- c("DESeq2", "t-test", "t-test")
# profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
# x=list(rnaseq, imputed_methyl, rppa)
# 
# # model_name -> res_pa_GMR_d_18_28.RData
# # id -> result18_28_GMR_d
# result_name <- paste(c('result',id,'_GMR_d'), collapse = '')
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, lim,
#                id = result_name, prob = p, Gamma = g, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)
# 
# model <- fit.classification(y=y, samples = samples, id = result_name, datapath = datapath, respath = respath,
#                             profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                             nFolds = 5, numTops=50, iter = 10)
# 
# 
# model_path <- paste(c('data/model/res_pa_GMR_d_',id,'_LOOCV.RData'), collapse = '')
# 
# name <- paste(c('res_pa_GMR_d_', id, '_LOOCV'), collapse='')
# assign(x = name, value = model)
# 
# save(list=name, file=file.path(model_path))

##############################################################
load('data/model/res_pa_GMR_d_18_17_LOOCV.RData.RData')
res_gmr_d_25 <- c(res_gmr_d_25, list(res_pa_GMR_d_18_17_LOOCV))

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
