# integrative DRW on combined feature data (updated in 2018/06/17)
# concat directed pathway graphs within each profile (GM & GMR & GMR_d & GMP)

# All gene symbols are converted to Entrez gene id

# Reduced PPI network are used
# Extract the protein(node) which is corresponding to RPPA protein and its neighbor node
# PPI relation : controls-phosphorylation-of

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)

#---------------------------------Result 10_2---------------------------------#
# PPI relation : controls-phosphorylation-of
dppipath <- file.path(datapath, 'DppiGraph_2.rda')
load(file.path(dppipath))

p <- DppiGraph_2
V(p)$name <-paste("p",V(p)$name,sep="")


#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_2", 
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_10_2 <- fit.classification(y=y, samples = samples, id="result10_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                     method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_10_2, file=file.path('data/model/res_pa_GM_10_2.RData'))

summary(res_pa_GM_10_2)
print(res_pa_GM_10_2$results)
print(res_pa_GM_10_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GM_10_2, id="result10_2", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_2", 
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_10_2 <- fit.classification(y=y, samples = samples, id="result10_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_10_2, file=file.path('data/model/res_pa_GMR_10_2.RData'))

summary(res_pa_GMR_10_2)
print(res_pa_GMR_10_2$results)
print(res_pa_GMR_10_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_10_2, id="result10_2", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_2", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_10_2 <- fit.classification(y=y, samples = samples, id="result10_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_10_2, file=file.path('data/model/res_pa_GMR_d_10_2.RData'))

summary(res_pa_GMR_d_10_2)
print(res_pa_GMR_d_10_2$results)
print(res_pa_GMR_d_10_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_d_10_2, id="result10_2", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_2", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_2 <- fit.classification(y=y, samples = samples, id="result10_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_2, file=file.path('data/model/res_pa_GMP_10_2.RData'))

summary(res_pa_GMP_10_2)
print(res_pa_GMP_10_2$results)
print(res_pa_GMP_10_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_2, id="result10_2", profile_name=profile_name, method="DRW", respath=respath)


# plot
title <- c("Result 10_2")
xlabs <- c("GM", "GMR", "GMR_d", "GMP")
res_models <- list(res_pa_GM_10_2, res_pa_GMR_10_2, res_pa_GMR_d_10_2, res_pa_GMP_10_2)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_10_2$results$Accuracy[1])