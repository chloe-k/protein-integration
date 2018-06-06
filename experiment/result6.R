# DRW-based single pathway profile on feature data (updated in 2018/05/31)
# directed pathway graphs was used in (G, M, R) model and PPI graph was used in (P) model

# Classifier : rf(Random Forest)

#------------------------- RNAseq(pathway) -------------------------#
testStatistic <- c("DESeq2")
profile_name <- c("rna")
x=list(rnaseq)

fit.iDRWPClass(x=x, y=y, globalGraph=g, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "G", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_G_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                               method = "DRW", pranking = "t-test", classifier = "rf",
                               nFolds = 5, numTops=50, iter = 50)

save(res_pa_G_RF, file=file.path('data/model/res_pa_G_RF.RData'))

summary(res_pa_G_RF)
print(res_pa_G_RF$results)
print(res_pa_G_RF$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_G_RF, profile_name=profile_name, method="DRW", respath=respath)

#------------------------- Methyl(pathway) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("meth")
x=list(imputed_methyl)


fit.iDRWPClass(x=x, y=y, globalGraph=m, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "M", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_M_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                               method = "DRW", pranking = "t-test", classifier = "rf",
                               nFolds = 5, numTops=50, iter = 50)

save(res_pa_M_RF, file=file.path('data/model/res_pa_M_RF.RData'))

summary(res_pa_M_RF)
print(res_pa_M_RF$results)
print(res_pa_M_RF$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_M_RF, profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using Pathway Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=r, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "R", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_R_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                               method = "DRW", pranking = "t-test", classifier = "rf",
                               nFolds = 5, numTops=50, iter = 50)

save(res_pa_R_RF, file=file.path('data/model/res_pa_R_RF.RData'))

summary(res_pa_R_RF)
print(res_pa_R_RF$results)
print(res_pa_R_RF$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_R_RF, profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using PPI Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=p, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "P", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_P_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                               method = "DRW", pranking = "t-test", classifier = "rf",
                               nFolds = 5, numTops=50, iter = 50)

save(res_pa_P_RF, file=file.path('data/model/res_pa_P_RF.RData'))

summary(res_pa_P_RF)
print(res_pa_P_RF$results)
print(res_pa_P_RF$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_P_RF, profile_name=profile_name, method="DRW", respath=respath)

# plot
xlabs <- c("G", "M", "R", "P")
res_models <- list(res_pa_G_RF, res_pa_M_RF, res_pa_R_RF, res_pa_P_RF)

perf_boxplot(xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_G_RF$results$Accuracy[1])