# In Result A, optimized hyperparameter was used for Random Walk with Restart(RWR) algorithm

# integrative DRW on combined feature data (updated in 2018/07/08)
# concat directed pathway graphs within each profile (GM & GMR & GMR_d & GMP)

# p=0.8 which is restart probability for network diffusion in diffus_ppi.R
# Gamma=  which is restart probability for DRW in DRW.R

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.

# DppiGraph(Entrez).rda was used

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result_A",
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_A <- fit.classification(y=y, samples = samples, id = "result_A", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_RF_A, file=file.path('data/model/res_pa_GM_A.RData'))

summary(res_pa_GM_A)
print(res_pa_GM_A$results)
print(res_pa_GM_A$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GM_A, id = "result_A", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result_A",
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_A <- fit.classification(y=y, samples = samples, id = "result_A", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_A, file=file.path('data/model/res_pa_GMR_A.RData'))

summary(res_pa_GMR_A)
print(res_pa_GMR_A$results)
print(res_pa_GMR_A$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_A, id = "result_A", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result_A",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_A <- fit.classification(y=y, samples = samples, id = "result_A", datapath = datapath, respath = respath, profile_name = profile_name,
                                         method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_A, file=file.path('data/model/res_pa_GMR_d_A.RData'))

summary(res_pa_GMR_d_A)
print(res_pa_GMR_d_A$results)
print(res_pa_GMR_d_A$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_d_A, id = "result_A", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result_A",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_A <- fit.classification(y=y, samples = samples, id = "result_A", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_A, file=file.path('data/model/res_pa_GMP_A.RData'))

summary(res_pa_GMP_A)
print(res_pa_GMP_A$results)
print(res_pa_GMP_A$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_A, id = "result_A", profile_name=profile_name, method="DRW", respath=respath)


# plot
title <- c("Result A")
xlabs <- c("GM_0.2", "GMR_0.2", "GMR_d_0.2", "GMP_0.2")
res_models <- list(res_pa_GM_20_0.7, res_pa_GMR_20_0.7, res_pa_GMR_d_20_0.7, res_pa_GMP_20_0.7)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, mean(res_pa_GM_20_0.7$results$Accuracy))
