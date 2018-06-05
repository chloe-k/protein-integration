# integrative DRW on combined feature data (updated in 2018/05/30)
# concat directed pathway graphs within each profile (GM & GMR & GMP)

# edge direction
# m -> g
# p -> g

# svmLinear classifier was used and remainder are same with result1.
# 'svmLinear' is Support Vector Machines with Linear Kernel.

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna", "meth")
x=list(rnaseq, imputed_methyl)

res_pa_GM_SVM <- fit.iDRWPClass(x=x, y=y, globalGraph=gm,
                            testStatistic= testStatistic, profile_name = profile_name,
                            datapath = datapath, respath = respath, pathSet=pathSet,
                            method = "DRW", samples = samples, pranking = "t-test", mode = "GM",
                            classifier = "svmLinear", nFolds = 5, iter = 50, AntiCorr=FALSE, DEBUG=TRUE)

save(res_pa_GM_SVM, file=file.path('data/model/res_pa_GM_SVM.RData'))

summary(res_pa_GM_SVM)
print(res_pa_GM_SVM$results)
print(res_pa_GM_SVM$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GM_SVM, profile_name=profile_name, method="DRW", classifier = "SVM", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa(Pathway_Graph)")
x=list(rnaseq, imputed_methyl, rppa)

res_pa_GMR_SVM <- fit.iDRWPClass(x=x, y=y, globalGraph=gmr,
                             testStatistic= testStatistic, profile_name = profile_name,
                             datapath = datapath, respath = respath, pathSet=pathSet,
                             method = "DRW", samples = samples, pranking = "t-test", mode = "GMR",
                             classifier = "svmLinear", nFolds = 5, iter = 50, AntiCorr=FALSE, DEBUG=TRUE)

save(res_pa_GMR_SVM, file=file.path('data/model/res_pa_GMR_SVM.RData'))

summary(res_pa_GMR_SVM)
print(res_pa_GMR_SVM$results)
print(res_pa_GMR_SVM$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_SVM, profile_name=profile_name, method="DRW", classifier = "SVM", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa(diffused_Pathway_Graph)")
x=list(rnaseq, imputed_methyl, rppa)

res_pa_GMR_d_SVM <- fit.iDRWPClass(x=x, y=y, globalGraph=gmr,
                                  testStatistic= testStatistic, profile_name = profile_name,
                                  datapath = datapath, respath = respath, pathSet=pathSet,
                                  method = "DRW", samples = samples, pranking = "t-test", mode = "GMP",
                                  classifier = 'svmLinear', nFolds = 5, iter = 50, AntiCorr=FALSE, DEBUG=TRUE)

save(res_pa_GMR_d_SVM, file=file.path('data/model/res_pa_GMR_d_SVM.RData'))

summary(res_pa_GMR_d_SVM)
print(res_pa_GMR_d_SVM$results)
print(res_pa_GMR_d_SVM$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_d_SVM, profile_name=profile_name, method="DRW", classifier = "SVM", respath=respath)



#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

res_pa_GMP_SVM <- fit.iDRWPClass(x=x, y=y, globalGraph=gmp,
                             testStatistic= testStatistic, profile_name = profile_name,
                             datapath = datapath, respath = respath, pathSet=pathSet,
                             method = "DRW", samples = samples, pranking = "t-test", mode = "GMP",
                             classifier = "svmLinear", nFolds = 5, iter = 50, AntiCorr=FALSE, DEBUG=TRUE)


save(res_pa_GMP_SVM, file=file.path('data/model/res_pa_GMP_SVM.RData'))

summary(res_pa_GMP_SVM)
print(res_pa_GMP_SVM$results)
print(res_pa_GMP_SVM$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_SVM, profile_name=profile_name, method="DRW", classifier = "SVM", respath=respath)


# plot
xlabs <- c("GM", "GMR", "GMP")
res_models <- list(res_pa_GM_SVM, res_pa_GMR_SVM, res_pa_GMR_d_SVM, res_pa_GMP_SVM)

perf_boxplot(xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_SVM$results$Accuracy)