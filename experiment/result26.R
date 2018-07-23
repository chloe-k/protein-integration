# integrative DRW on combined feature data (updated in 2018/07/20)
# concat directed pathway graphs within each profile (GMR)

# All gene symbols are converted to Entrez gene id
# 5-fold CV(10 iters) was performed for tuning parameter in Random Forest.
# 5-fold CV(10 iters) was performed for get top N pathways.
# LOOCV was performed for model evaluation

# Dppigraph(Entrez).rda was used 
# Dup_rppa_data.RData was used in GM, GMR -> unmapped genes with rppa proteins are removed in rnaseq, imputed_methyl profiles
# Entrez_data.RData was used in GM_base

# Gamma = 0.8 was used

# data dimension
# rnaseq - 183 genes, 376 samples
# imputed_methyl - 176 genes, 376 samples
# rppa - 184 genes, 376 samples

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
               id = "result26_GM", prob = 0.001, Gamma = 0.8, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_26 <- fit.classification(y=y, samples = samples, id = "result26_GM", datapath = datapath, respath = respath,
                                   profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                   nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_26, file=file.path('data/model/res_pa_GM_26.RData'))

write.SigFeatures(res_fit=res_pa_GM_26, id="result26_GM", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_26))

#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(res26)", "meth(res26)", "rppa(res26)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result26_GMR", prob = 0.001, Gamma = 0.8, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_26 <- fit.classification(y=y, samples = samples, id = "result26_GMR", datapath = datapath, respath = respath,
                                    profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                    nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_26, file=file.path('data/model/res_pa_GMR_26.RData'))

write.SigFeatures(res_fit=res_pa_GMR_26, id="result26_GMR", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_26))

#------------------------- RNAseq + Methyl (baseline) -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result26_base_GM", prob = 0.001, Gamma = 0.8, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_26_base <- fit.classification(y=y, samples = samples, id = "result26_base_GM", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_26_base, file=file.path('data/model/res_pa_GM_26_base.RData'))

write.SigFeatures(res_fit=res_pa_GM_26_base, id="result26_base_GM", profile_name=profile_name, method="DRW", respath=respath)

baseline <- max(res_pa_GM_26_base$results$Accuracy)

# Plot for Result26

title <- c("Result 26")
xlabs <- c("GM", "GMR")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.05, perf_max = perf_max+0.05, baseline)
