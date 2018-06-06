# integrative DRW on combined feature data (updated in 2018/06/02)
# concat directed pathway graphs within each profile (GM & GMR & GMR_d & GMP)

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna", "meth")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                method = "DRW", pranking = "t-test", classifier = "rf",
                                nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_RF, file=file.path('data/model/res_pa_GM_RF.RData'))

summary(res_pa_GM_RF)
print(res_pa_GM_RF$results)
print(res_pa_GM_RF$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GM_RF, profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa(Pathway_Graph)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                 method = "DRW", pranking = "t-test", classifier = "rf",
                                 nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_RF, file=file.path('data/model/res_pa_GMR_RF.RData'))

summary(res_pa_GMR_RF)
print(res_pa_GMR_RF$results)
print(res_pa_GMR_RF$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_RF, profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa(diffused_Pathway_Graph)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_RF, file=file.path('data/model/res_pa_GMR_d_RF.RData'))

summary(res_pa_GMR_d_RF)
print(res_pa_GMR_d_RF$results)
print(res_pa_GMR_d_RF$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_d_RF, profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                 method = "DRW", pranking = "t-test", classifier = "rf",
                                 nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_RF, file=file.path('data/model/res_pa_GMP_RF.RData'))

summary(res_pa_GMP_RF)
print(res_pa_GMP_RF$results)
print(res_pa_GMP_RF$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_RF, profile_name=profile_name, method="DRW", respath=respath)


# plot
xlabs <- c("GM", "GMR", "GMR_d", "GMP")
res_models <- list(res_pa_GM_RF, res_pa_GMR_RF, res_pa_GMR_d_RF, res_pa_GMP_RF)

perf_boxplot(xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_RF$results$Accuracy[1])
