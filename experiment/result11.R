# DRW-based single pathway profile on feature data (updated in 2018/06/25)
# directed pathway graphs was used in (G, M, R) model and PPI graph was used in (P) model

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.

# Dppigraph(Entrez).rda was used

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)

#------------------------- RNAseq(pathway) -------------------------#
testStatistic <- c("DESeq2")
profile_name <- c("rna(Entrez)")
x=list(rnaseq)

fit.iDRWPClass(x=x, y=y, globalGraph=g, testStatistic= testStatistic, profile_name = profile_name, id = "result11_G",
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "G", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_G_11 <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "DRW", pranking = "t-test", classifier = "rf", id = "result11_G",
                                  nFolds = 5, numTops=50, iter = 50)

save(res_pa_G_11, file=file.path('data/model/res_pa_G_11.RData'))

write.SigFeatures(res_fit=res_pa_G_11, id = "result11_G", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- Methyl(pathway) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("meth(Entrez)")
x=list(imputed_methyl)


fit.iDRWPClass(x=x, y=y, globalGraph=m, testStatistic= testStatistic, profile_name = profile_name, id = "result11_M",
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "M", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_M_11 <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "DRW", pranking = "t-test", classifier = "rf", id = "result11_M",
                                  nFolds = 5, numTops=50, iter = 50)

save(res_pa_M_11, file=file.path('data/model/res_pa_M_11.RData'))

write.SigFeatures(res_fit=res_pa_M_11, id = "result11_M", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using Pathway Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa(Entrez)")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=r, testStatistic= testStatistic, profile_name = profile_name, id = "result11_R",
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "R", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_R_11 <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "DRW", pranking = "t-test", classifier = "rf", id = "result11_R",
                                  nFolds = 5, numTops=50, iter = 50)

save(res_pa_R_11, file=file.path('data/model/res_pa_R_11.RData'))

write.SigFeatures(res_fit=res_pa_R_11, id = "result11_R", profile_name=profile_name, method="DRW", respath=respath)

# #------------------------- RPPA(pathway using PPI Graph) -------------------------#
# testStatistic <- c("t-test")
# profile_name <- c("rppa(Entrez)")
# x=list(rppa)
# 
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=p, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
#                pranking = "t-test", mode = "P", AntiCorr=FALSE, DEBUG=TRUE)
# 
# res_pa_P_RF_11 <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
#                                   method = "DRW", pranking = "t-test", classifier = "rf",
#                                   nFolds = 5, numTops=50, iter = 50)
# 
# save(res_pa_P_RF_11, file=file.path('data/model/res_pa_P_RF_11.RData'))
# 
# summary(res_pa_P_RF_11)
# print(res_pa_P_RF_11$results)
# print(res_pa_P_RF_11$resample$Accuracy)
# 
# write.SigFeatures(res_fit=res_pa_P_RF_11, profile_name=profile_name, method="DRW", respath=respath)

# plot
title <- c("Result 11")
xlabs <- c("G", "M", "R")
res_models <- list(res_pa_G_RF_11, res_pa_M_RF_11, res_pa_R_RF_11, res_pa_P_RF_11)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_G_RF_11$results$Accuracy[1])