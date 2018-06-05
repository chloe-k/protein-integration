# integrative DRW on combined feature data (updated in 2018/05/30)
# concat directed pathway graphs within each profile (GM & GMR & GMP)

# edge direction
# m -> g
# p -> g

# Classifier : glm(Generalized Linear Model)

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna", "meth")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                                 method = "DRW", pranking = "t-test", classifier = "glm",
                                 nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM, file=file.path('data/model/res_pa_GM.RData'))

summary(res_pa_GM)
print(res_pa_GM$results)
print(res_pa_GM$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GM, profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa(Pathway_Graph)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                                 method = "DRW", pranking = "t-test", classifier = "glm",
                                 nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR, file=file.path('data/model/res_pa_GMR.RData'))

summary(res_pa_GMR)
print(res_pa_GMR$results)
print(res_pa_GMR$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR, profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa(diffused_Pathway_Graph)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                                   method = "DRW", pranking = "t-test", classifier = "glm",
                                   nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d, file=file.path('data/model/res_pa_GMR_d.RData'))

summary(res_pa_GMR_d)
print(res_pa_GMR_d$results)
print(res_pa_GMR_d$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_d, profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                                   method = "DRW", pranking = "t-test", classifier = "glm",
                                   nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP, file=file.path('data/model/res_pa_GMP.RData'))

summary(res_pa_GMP)
print(res_pa_GMP$results)
print(res_pa_GMP$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP, profile_name=profile_name, method="DRW", respath=respath)


# plot
xlabs <- c("GM", "GMR", "GMR_d", "GMP")
res_models <- list(res_pa_GM, res_pa_GMR, res_pa_GMR_d, res_pa_GMP)

perf_boxplot(xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM$results$Accuracy)