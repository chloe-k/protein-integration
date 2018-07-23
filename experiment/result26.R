# integrative DRW on combined feature data (updated in 2018/07/20)
# concat directed pathway graphs within each profile (GMR)

# All gene symbols are converted to Entrez gene id
# 5-fold CV(10 iters) was performed for tuning parameter in Random Forest.
# 5-fold CV(10 iters) was performed for get top N pathways.
# LOOCV was performed for model evaluation

# Dppigraph(Entrez).rda was used
# Dup_rppa_data.RData was used -> unmapped genes with rppa proteins are removed in rnaseq, imputed_methyl profiles

# Gamma = 0.4 was used

# data dimension
# rnaseq - 183 genes, 402 samples
# imputed_methyl - 176 genes, 402 samples
# rppa - 184 genes, 402 samples

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)




########################################### Result26 #######################################
res_models <- list()
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(res26)", "meth(res26)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result26_GM", prob = 0.001, Gamma = 0.4, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_26 <- fit.classification(y=y, samples = samples, id = "result26_GM", datapath = datapath, respath = respath,
                                   profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                   nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_26, file=file.path('data/model/res_pa_GM_26.RData'))

res_models <- c(res_models, list(res_pa_GM_26))

#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(res26)", "meth(res26)", "rppa(res26)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result26_GMR", prob = 0.001, Gamma = 0.4, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_26 <- fit.classification(y=y, samples = samples, id = "result26_GMR", datapath = datapath, respath = respath,
                                    profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                    nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_26, file=file.path('data/model/res_pa_GMR_26.RData'))

res_models <- c(res_models, list(res_pa_GMR_26))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <-  c("rna(res26)", "meth(res26)", "rppa(res26)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result26_GMR_d", prob = 0.4, Gamma = 0.4, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_26 <- fit.classification(y=y, samples = samples, id = "result26_GMR_d", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_26, file=file.path('data/model/res_pa_GMR_d_26.RData'))

res_models <- c(res_models, list(res_pa_GMR_d_26))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <-  c("rna(res26)", "meth(res26)", "rppa(res26)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result26_1_GMR_d", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_26_1 <- fit.classification(y=y, samples = samples, id = "result26_1_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_26_1, file=file.path('data/model/res_pa_GMR_d_26_1.RData'))

res_models <- c(res_models, list(res_pa_GMR_d_26_1))

############################################################################
# Plot for Result26

title <- c("Result 26")
xlabs <- c("GM", "GMR", "GMR_d", "GMR_d_1")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.05, perf_max = perf_max+0.05)