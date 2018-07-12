# integrative DRW on combined feature data (updated in 2018/07/08)
# concat directed pathway graphs within each profile (GM & GMR & GMR_d & GMP)

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.
# weighted PPI was constructed by STRING database
# combined_score/100 was applied as edge weight

# Dppigraph_W_str.rda was used
# restart probability(in diffus_ppi.R) is 0.8

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)

# #------------------------- RNAseq + Methyl -------------------------#
# gm <- g %du% m
# testStatistic <- c("DESeq2", "t-test")
# profile_name <- c("rna(Entrez)", "meth(Entrez)")
# x=list(rnaseq, imputed_methyl)
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result22",
#                pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)
# 
# res_pa_GM_RF_22 <- fit.classification(y=y, samples = samples, id = "result22", datapath = datapath, respath = respath, profile_name = profile_name,
#                                       method = "DRW", pranking = "t-test", classifier = "rf",
#                                       nFolds = 5, numTops=50, iter = 50)
# 
# save(res_pa_GM_RF_22, file=file.path('data/model/res_pa_GM_RF_22.RData'))
# 
# summary(res_pa_GM_RF_22)
# print(res_pa_GM_RF_22$results)
# 
# write.SigFeatures(res_fit=res_pa_GM_RF_22, id = "result22", profile_name=profile_name, method="DRW", respath=respath)
# 
# 
# #------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
# gmr <- g %du% m %du% r
# testStatistic <- c("DESeq2", "t-test", "t-test")
# profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
# x=list(rnaseq, imputed_methyl, rppa)
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result22",
#                pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)
# 
# res_pa_GMR_RF_22 <- fit.classification(y=y, samples = samples, id = "result22", datapath = datapath, respath = respath, profile_name = profile_name,
#                                        method = "DRW", pranking = "t-test", classifier = "rf",
#                                        nFolds = 5, numTops=50, iter = 50)
# 
# save(res_pa_GMR_RF_22, file=file.path('data/model/res_pa_GMR_RF_22.RData'))
# 
# summary(res_pa_GMR_RF_22)
# print(res_pa_GMR_RF_22$results)
# 
# write.SigFeatures(res_fit=res_pa_GMR_RF_22, id = "result22", profile_name=profile_name, method="DRW", respath=respath)
# 
# 
# #------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
# gmr <- list(g, m, r)
# testStatistic <- c("DESeq2", "t-test", "t-test")
# profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
# x=list(rnaseq, imputed_methyl, rppa)
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result22",
#                pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)
# 
# res_pa_GMR_d_RF_22 <- fit.classification(y=y, samples = samples, id = "result22", datapath = datapath, respath = respath, profile_name = profile_name,
#                                          method = "DRW", pranking = "t-test", classifier = "rf",
#                                          nFolds = 5, numTops=50, iter = 50)
# 
# save(res_pa_GMR_d_RF_22, file=file.path('data/model/res_pa_GMR_d_RF_22.RData'))
# 
# summary(res_pa_GMR_d_RF_22)
# print(res_pa_GMR_d_RF_22$results)
# 
# write.SigFeatures(res_fit=res_pa_GMR_d_RF_22, id = "result22", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
dppipath <- file.path(datapath, 'DppiGraph_W_str.rda')  # directed edge PPI
load(file.path(dppipath))

p <- DppiGraph
V(p)$name <-paste("p",V(p)$name,sep="")

testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result22",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_22 <- fit.classification(y=y, samples = samples, id = "result22", datapath = datapath, respath = respath, profile_name = profile_name,
                                    method = "DRW", pranking = "t-test", classifier = "rf",
                                    nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_22, file=file.path('data/model/res_pa_GMP_22.RData'))

summary(res_pa_GMP_22)
print(res_pa_GMP_22$results)

write.SigFeatures(res_fit=res_pa_GMP_22, id = "result22", profile_name=profile_name, method="DRW", respath=respath)


# plot
title <- c("Result 22")
xlabs <- c("GM", "GMR", "GMR_d", "GMP(unweighted STRING PPI)", "GMP(weighted STRING PPI)")
res_models <- list(res_pa_GM_RF_13, res_pa_GMR_RF_13, res_pa_GMR_d_RF_13, res_pa_GMP_21, res_pa_GMP_22)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_RF_13$results$Accuracy[1])
