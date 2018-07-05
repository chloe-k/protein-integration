# DRW-based single pathway profile on feature data (updated in 2018/06/25)
# directed pathway graphs was used in (G, M, R) model and PPI graph was used in (P) model

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.

# DppiGraph_rdc.rda was used

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)

#------------------------- RNAseq(pathway) -------------------------#
testStatistic <- c("DESeq2")
profile_name <- c("rna(Entrez)")
x=list(rnaseq)

fit.iDRWPClass(x=x, y=y, globalGraph=g, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result14", 
               pranking = "t-test", mode = "G", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_G_14 <- fit.classification(y=y, samples = samples, id="result14", datapath = datapath, respath = respath, profile_name = profile_name,
                                     method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_G_14, file=file.path('data/model/res_pa_G_14.RData'))

summary(res_pa_G_14)
print(res_pa_G_14$results)
print(res_pa_G_14$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_G_14, id="result14", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- Methyl(pathway) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("meth(Entrez)")
x=list(imputed_methyl)


fit.iDRWPClass(x=x, y=y, globalGraph=m, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result14", 
               pranking = "t-test", mode = "M", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_M_14 <- fit.classification(y=y, samples = samples, id="result14", datapath = datapath, respath = respath, profile_name = profile_name,
                                     method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_M_14, file=file.path('data/model/res_pa_M_14.RData'))

summary(res_pa_M_14)
print(res_pa_M_14$results)
print(res_pa_M_14$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_M_14, id="result14", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using Pathway Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa(Entrez)")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=r, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result14", 
               pranking = "t-test", mode = "R", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_R_14 <- fit.classification(y=y, samples = samples, id="result14", datapath = datapath, respath = respath, profile_name = profile_name,
                                     method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_R_14, file=file.path('data/model/res_pa_R_14.RData'))

summary(res_pa_R_14)
print(res_pa_R_14$results)
print(res_pa_R_14$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_R_14, id="result14", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using PPI Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa(Entrez)")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=p, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result14", 
               pranking = "t-test", mode = "P", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_P_14 <- fit.classification(y=y, samples = samples, id="result14", datapath = datapath, respath = respath, profile_name = profile_name,
                                     method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_P_14, file=file.path('data/model/res_pa_P_14.RData'))

summary(res_pa_P_14)
print(res_pa_P_14$results)
print(res_pa_P_14$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_P_14, id="result14", profile_name=profile_name, method="DRW", respath=respath)

# plot
title <- c("Result 14")
xlabs <- c("G", "M", "R", "P")
res_models <- list(res_pa_G_14, res_pa_M_14, res_pa_R_14, res_pa_P_14)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_G_14$results$Accuracy[1])