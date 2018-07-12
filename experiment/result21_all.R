# integrative DRW on combined feature data (updated in 2018/07/09)
# concat directed pathway graphs within each profile (GM & GMR & GMR_d & GMP)

# All gene symbols are converted to Entrez gene id

# Dppigraph_W_str.rda was used for GMP
# combined_score/1000 was applied as edge weight(Edge weight was only applied in model using P)

# restart probability(in diffus_ppi.R) p is 0.8
# restart probability(in DRW.R) Gamma is 0.1(In res_pa_*_A, Gamma is also 0.1)


# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)

#---------------------------------Result 21_1---------------------------------#
# PPI relation : interacts-with
dppipath <- file.path(datapath, 'DppiGraph21_1.rda')
load(file.path(dppipath))

p <- DppiGraph21_1
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result21_1",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_21_1 <- fit.classification(y=y, samples = samples, id="result21_1", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_21_1, file=file.path('data/model/res_pa_GMP_21_1.RData'))

summary(res_pa_GMP_21_1)
print(res_pa_GMP_21_1$results)

write.SigFeatures(res_fit=res_pa_GMP_21_1, id="result21_1", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 21_2---------------------------------#
# PPI relation : controls-phosphorylation-of
dppipath <- file.path(datapath, 'DppiGraph21_2.rda')
load(file.path(dppipath))

p <- DppiGraph21_2
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result21_2",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_21_2 <- fit.classification(y=y, samples = samples, id="result21_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_21_2, file=file.path('data/model/res_pa_GMP_21_2.RData'))

summary(res_pa_GMP_21_2)
print(res_pa_GMP_21_2$results)

write.SigFeatures(res_fit=res_pa_GMP_21_2, id="result21_2", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 21_3---------------------------------#
dppipath <- file.path(datapath, 'DppiGraph21_3.rda')
load(file.path(dppipath))

p <- DppiGraph21_3
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result21_3",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_21_3 <- fit.classification(y=y, samples = samples, id="result21_3", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_21_3, file=file.path('data/model/res_pa_GMP_21_3.RData'))

summary(res_pa_GMP_21_3)
print(res_pa_GMP_21_3$results)

write.SigFeatures(res_fit=res_pa_GMP_21_3, id="result21_3", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 21_4---------------------------------#

dppipath <- file.path(datapath, 'DppiGraph21_4.rda')
load(file.path(dppipath))

p <- DppiGraph21_4
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result21_4",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_21_4 <- fit.classification(y=y, samples = samples, id="result21_4", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_21_4, file=file.path('data/model/res_pa_GMP_21_4.RData'))

summary(res_pa_GMP_21_4)
print(res_pa_GMP_21_4$results)

write.SigFeatures(res_fit=res_pa_GMP_21_4, id="result21_4", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 21_5---------------------------------#

dppipath <- file.path(datapath, 'DppiGraph21_5.rda')
load(file.path(dppipath))

p <- DppiGraph21_5
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result21_5",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_21_5 <- fit.classification(y=y, samples = samples, id="result21_5", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_21_5, file=file.path('data/model/res_pa_GMP_21_5.RData'))

summary(res_pa_GMP_21_5)
print(res_pa_GMP_21_5$results)

write.SigFeatures(res_fit=res_pa_GMP_21_5, id="result21_5", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 21_6---------------------------------#

dppipath <- file.path(datapath, 'DppiGraph21_6.rda')
load(file.path(dppipath))

p <- DppiGraph21_6
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result21_6",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_21_6 <- fit.classification(y=y, samples = samples, id="result21_6", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_21_6, file=file.path('data/model/res_pa_GMP_21_6.RData'))

summary(res_pa_GMP_21_6)
print(res_pa_GMP_21_6$results)

write.SigFeatures(res_fit=res_pa_GMP_21_6, id="result21_6", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 21_7---------------------------------#

dppipath <- file.path(datapath, 'DppiGraph21_7.rda')
load(file.path(dppipath))

p <- DppiGraph21_7
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result21_7",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_21_7 <- fit.classification(y=y, samples = samples, id="result21_7", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_21_7, file=file.path('data/model/res_pa_GMP_21_7.RData'))

summary(res_pa_GMP_21_7)
print(res_pa_GMP_21_7$results)

write.SigFeatures(res_fit=res_pa_GMP_21_7, id="result21_7", profile_name=profile_name, method="DRW", respath=respath)



# plot
title <- c("Result 21_all")
xlabs <- c("GM", "GMR", "GMR_d", "GMP", "GMP21_1", "GMP21_2", "GMP21_3", "GMP21_4", "GMP21_5", "GMP21_6", "GMP21_7")
res_models <- list(res_pa_GM_A, res_pa_GMR_A, res_pa_GMR_d_A, res_pa_GMP_A, res_pa_GMP_21_1, res_pa_GMP_21_2, res_pa_GMP_21_3, res_pa_GMP_21_4, res_pa_GMP_21_5, res_pa_GMP_21_6, res_pa_GMP_21_7)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, mean(res_pa_GM_A$results$Accuracy))