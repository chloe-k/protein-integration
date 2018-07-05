# integrative DRW on combined feature data (updated in 2018/06/17)
# concat directed pathway graphs within each profile (GM & GMR & GMR_d & GMP)

# All gene symbols are converted to Entrez gene id

# Reduced PPI network are used
# Extract the protein(node) which is corresponding to RPPA protein and its neighbor node
# PPI relation : interacts-with

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)

#---------------------------------Result 10_1---------------------------------#
# PPI relation : interacts-with
dppipath <- file.path(datapath, 'DppiGraph_1.rda')
load(file.path(dppipath))

p <- DppiGraph_1
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_1", 
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_10_1 <- fit.classification(y=y, samples = samples, id="result10_1", datapath = datapath, respath = respath, profile_name = profile_name,
                                     method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_10_1, file=file.path('data/model/res_pa_GM_10_1.RData'))

summary(res_pa_GM_10_1)
print(res_pa_GM_10_1$results)
print(res_pa_GM_10_1$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GM_10_1, id="result10_1", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_1", 
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_10_1 <- fit.classification(y=y, samples = samples, id="result10_1", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_10_1, file=file.path('data/model/res_pa_GMR_10_1.RData'))

summary(res_pa_GMR_10_1)
print(res_pa_GMR_10_1$results)
print(res_pa_GMR_10_1$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_10_1, id="result10_1", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_1", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_10_1 <- fit.classification(y=y, samples = samples, id="result10_1", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_10_1, file=file.path('data/model/res_pa_GMR_d_10_1.RData'))

summary(res_pa_GMR_d_10_1)
print(res_pa_GMR_d_10_1$results)
print(res_pa_GMR_d_10_1$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_d_10_1, id="result10_1", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_1", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_1 <- fit.classification(y=y, samples = samples, id="result10_1", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_1, file=file.path('data/model/res_pa_GMP_10_1.RData'))

summary(res_pa_GMP_10_1)
print(res_pa_GMP_10_1$results)
print(res_pa_GMP_10_1$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_1, id="result10_1", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_2---------------------------------#
# PPI relation : controls-phosphorylation-of
dppipath <- file.path(datapath, 'DppiGraph_2.rda')
load(file.path(dppipath))

p <- DppiGraph_2
V(p)$name <-paste("p",V(p)$name,sep="")

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

#---------------------------------Result 10_3---------------------------------#
# PPI relation : catalysis-precedes
dppipath <- file.path(datapath, 'DppiGraph_3.rda')
load(file.path(dppipath))

p <- DppiGraph_3
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_3", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_3 <- fit.classification(y=y, samples = samples, id="result10_3", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_3, file=file.path('data/model/res_pa_GMP_10_3.RData'))

summary(res_pa_GMP_10_3)
print(res_pa_GMP_10_3$results)
print(res_pa_GMP_10_3$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_3, id="result10_3", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_4---------------------------------#
# PPI relation : controls-expression-of
dppipath <- file.path(datapath, 'DppiGraph_4.rda')
load(file.path(dppipath))

p <- DppiGraph_4
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_4", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_4 <- fit.classification(y=y, samples = samples, id="result10_4", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_4, file=file.path('data/model/res_pa_GMP_10_4.RData'))

summary(res_pa_GMP_10_4)
print(res_pa_GMP_10_4$results)
print(res_pa_GMP_10_4$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_4, id="result10_4", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_5---------------------------------#
# PPI relation : controls-state-change-of
dppipath <- file.path(datapath, 'DppiGraph_5.rda')
load(file.path(dppipath))

p <- DppiGraph_5
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_5", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_5 <- fit.classification(y=y, samples = samples, id="result10_5", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_5, file=file.path('data/model/res_pa_GMP_10_5.RData'))

summary(res_pa_GMP_10_5)
print(res_pa_GMP_10_5$results)
print(res_pa_GMP_10_5$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_5, id="result10_5", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_6---------------------------------#
# PPI relation : in-complex-with
dppipath <- file.path(datapath, 'DppiGraph_6.rda')
load(file.path(dppipath))

p <- DppiGraph_6
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_6", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_6 <- fit.classification(y=y, samples = samples, id="result10_6", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_6, file=file.path('data/model/res_pa_GMP_10_6.RData'))

summary(res_pa_GMP_10_6)
print(res_pa_GMP_10_6$results)
print(res_pa_GMP_10_6$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_6, id="result10_6", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_7---------------------------------#
# PPI relation : controls-transport-of
dppipath <- file.path(datapath, 'DppiGraph_7.rda')
load(file.path(dppipath))

p <- DppiGraph_7
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_7", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_7 <- fit.classification(y=y, samples = samples, id="result10_7", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_7, file=file.path('data/model/res_pa_GMP_10_7.RData'))

summary(res_pa_GMP_10_7)
print(res_pa_GMP_10_7$results)
print(res_pa_GMP_10_7$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_7, id="result10_7", profile_name=profile_name, method="DRW", respath=respath)



# plot
title <- c("Result 10_all")
xlabs <- c("GM", "GMR", "GMR_d", "GMP_1", "GMP_2", "GMP_3", "GMP_4", "GMP_5", "GMP_6", "GMP_7")
res_models <- list(res_pa_GM_10_1, res_pa_GMR_10_1, res_pa_GMR_d_10_1, res_pa_GMP_10_1, res_pa_GMP_10_2, res_pa_GMP_10_3, res_pa_GMP_10_4, res_pa_GMP_10_5, res_pa_GMP_10_6, res_pa_GMP_10_7)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_10_1$results$Accuracy[1])