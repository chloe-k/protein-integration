# integrative DRW on combined feature data (updated in 2018/07/10)
# concat directed pathway graphs within each profile (GM & GMR & GMR_d & GMP)

# For PPI network diffusion, Random Walk with Restart(RWR) algorithm was used.
# In order to find optimized restart probability in PPI diffusion.
# Grid search was performed about combination of p=[0.001, 0.01, 0.2, 0.4, 0.6, 0.8] and Gamma=[0, 0.2, 0.4, 0.6, 0.8]
# p=0.5 had used in before

# parameter tuning for GM model, extra experiment was performed by adding Gamma = [0.7, 0.75, 0.85, 0.9, 0.95]

# All gene symbols are converted to Entrez gene id
# 5-fold CV(20 iters) was performed.

# Dppigraph(Entrez).rda was used

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)



################################################### Result24 #################################################

res_models <- list()
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result24_GM", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_24 <- fit.classification(y=y, samples = samples, id = "result24_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     numTops=50)

save(res_pa_GM_24, file=file.path('data/model/res_pa_GM_24.RData'))

res_models <- c(res_models, list(res_pa_GM_24))

#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result24_GMR", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_24 <- fit.classification(y=y, samples = samples, id = "result24_GMR", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      numTops=50)

save(res_pa_GMR_24, file=file.path('data/model/res_pa_GMR_24.RData'))

res_models <- c(res_models, list(res_pa_GMR_24))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result24_GMP", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_24 <- fit.classification(y=y, samples = samples, id = "result24_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      numTops=50)


save(res_pa_GMP_24, file=file.path('data/model/res_pa_GMP_24.RData'))

res_models <- c(res_models, list(res_pa_GMP_24))


#########################################################################################################################################
# plot


# Plot for GM models
title <- c("Result 24 - all graphs are diffused")
xlabs <- c("GM_d", "GMR_d", "GMP_d")

perf_min <- min(sapply(X = res_models, FUN = function(x){mean(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){mean(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.02, perf_max = perf_max+0.02)