# DRW-based single pathway profile on feature data (updated in 2018/05/31)
# directed pathway graphs was used in (G, M, R) model and PPI graph was used in (P) model

# Classifier : glm(Generalized Linear Model)

#------------------------- RNAseq(pathway) -------------------------#
testStatistic <- c("DESeq2")
profile_name <- c("rna")
x=list(rnaseq)

fit.iDRWPClass(x=x, y=y, globalGraph=g, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "G", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_G <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                   method = "DRW", pranking = "t-test", classifier = "glm",
                                   nFolds = 5, numTops=50, iter = 50)

save(res_pa_G, file=file.path('data/model/res_pa_G.RData'))

summary(res_pa_G)
print(res_pa_G$results)
print(res_pa_G$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_G, profile_name=profile_name, method="DRW", respath=respath)

#------------------------- Methyl(pathway) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("meth")
x=list(imputed_methyl)


fit.iDRWPClass(x=x, y=y, globalGraph=m, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "M", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_M <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                               method = "DRW", pranking = "t-test", classifier = "glm",
                               nFolds = 5, numTops=50, iter = 50)

save(res_pa_M, file=file.path('data/model/res_pa_M.RData'))

summary(res_pa_M)
print(res_pa_M$results)
print(res_pa_M$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_M, profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using Pathway Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=r, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "R", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_R <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                               method = "DRW", pranking = "t-test", classifier = "glm",
                               nFolds = 5, numTops=50, iter = 50)

save(res_pa_R, file=file.path('data/model/res_pa_R.RData'))

summary(res_pa_R)
print(res_pa_R$results)
print(res_pa_R$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_R, profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using PPI Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=p, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "P", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_P <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                               method = "DRW", pranking = "t-test", classifier = "glm",
                               nFolds = 5, numTops=50, iter = 50)

save(res_pa_P, file=file.path('data/model/res_pa_P.RData'))

summary(res_pa_P)
print(res_pa_P$results)
print(res_pa_P$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_P, profile_name=profile_name, method="DRW", respath=respath)

# plot
xlabs <- c("G", "M", "R", "P")
res_models <- list(res_pa_G, res_pa_M, res_pa_R, res_pa_P)

perf_boxplot(xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_G$results$Accuracy)