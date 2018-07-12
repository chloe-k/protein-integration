# DRW-based single pathway profile on feature data (updated in 2018/07/08)
# directed pathway graphs was used as (G, M, R) model and PPI graph(Pathway Commons) was used as (P) model

# All KEGG graph/PPI graph was diffused by corresponding data
# restart probability(in diffus_ppi.R) is 0.8

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.

# DppiGraph(Entrez).rda was used

# Classifier : rf(Random Forest)

#------------------------- RNAseq(pathway) -------------------------#
testStatistic <- c("DESeq2")
profile_name <- c("rna(Entrez)")
x=list(rnaseq)

fit.iDRWPClass(x=x, y=y, globalGraph=g, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result18_2", 
               pranking = "t-test", mode = "G", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_G_18_2 <- fit.classification(y=y, samples = samples, id="result18_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "DRW", pranking = "t-test", classifier = "rf",
                                  nFolds = 5, numTops=50, iter = 50)

save(res_pa_G_18_2, file=file.path('data/model/res_pa_G_18_2.RData'))

summary(res_pa_G_18_2)
print(res_pa_G_18_2$results)
print(res_pa_G_18_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_G_18_2, id="result18_2", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- Methyl(pathway) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("meth(Entrez)")
x=list(imputed_methyl)


fit.iDRWPClass(x=x, y=y, globalGraph=m, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result18_2", 
               pranking = "t-test", mode = "M", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_M_18_2 <- fit.classification(y=y, samples = samples, id="result18_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "DRW", pranking = "t-test", classifier = "rf",
                                  nFolds = 5, numTops=50, iter = 50)

save(res_pa_M_18_2, file=file.path('data/model/res_pa_M_18_2.RData'))

summary(res_pa_M_18_2)
print(res_pa_M_18_2$results)
print(res_pa_M_18_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_M_18_2, id="result18_2", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using Pathway Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa(Entrez)")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=r, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result18_2", 
               pranking = "t-test", mode = "R", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_R_18_2 <- fit.classification(y=y, samples = samples, id="result18_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "DRW", pranking = "t-test", classifier = "rf",
                                  nFolds = 5, numTops=50, iter = 50)

save(res_pa_R_18_2, file=file.path('data/model/res_pa_R_18_2.RData'))

summary(res_pa_R_18_2)
print(res_pa_R_18_2$results)
print(res_pa_R_18_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_R_18_2, id="result18_2", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using PPI Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa(Entrez)")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=p, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result18_2", 
               pranking = "t-test", mode = "P", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_P_18_2 <- fit.classification(y=y, samples = samples, id="result18_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "DRW", pranking = "t-test", classifier = "rf",
                                  nFolds = 5, numTops=50, iter = 50)

save(res_pa_P_18_2, file=file.path('data/model/res_pa_P_18_2.RData'))

summary(res_pa_P_18_2)
print(res_pa_P_18_2$results)
print(res_pa_P_18_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_P_18_2, id="result18_2", profile_name=profile_name, method="DRW", respath=respath)

# plot
title <- c("Result 18_2")
xlabs <- c("G", "M", "R", "P")
res_models <- list(res_pa_G_18_2, res_pa_M_18_2, res_pa_R_18_2, res_pa_P_18_2)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_G_18_2$results$Accuracy[1])