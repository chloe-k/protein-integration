# integrative DRW on combined feature data (updated in 2018/07/07)
# concat directed pathway graphs within each profile (GM & GMR & GMR_d & GMP)

# For PPI network diffusion, Heat Diffusion Process on a Laplacian Matrix was performed.
# In order to find optimized time point when heat is measured in PPI diffusion.
# Grid search was performed in [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9]


# All gene symbols are converted to Entrez gene id
# LOOCV was performed.

# Dppigraph(Entrez).rda was used

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)


res_models <- list()
load('data/model/res_pa_GM_RF_13.RData')
load('data/model/res_pa_GMR_RF_13.RData')
res_models <- c(res_models, list(res_pa_GM_RF_13, res_pa_GMR_RF_13))

###########################################################  Prob = 0.01  #######################################################
prob <- 0.01
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.01 <- fit.classification(y=y, samples = samples, id = "result19", datapath = datapath, respath = respath, profile_name = profile_name,
                                           method = "DRW", pranking = "t-test", classifier = "rf",
                                           nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.01, file=file.path('data/model/res_pa_GMR_d_19_0.01.RData'))

summary(res_pa_GMR_d_19_0.01)
print(res_pa_GMR_d_19_0.01$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.01, id = "result19", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.01 <- fit.classification(y=y, samples = samples, id = "result19", datapath = datapath, respath = respath, profile_name = profile_name,
                                         method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.01, file=file.path('data/model/res_pa_GMP_19_0.01.RData'))

summary(res_pa_GMP_19_0.01)
print(res_pa_GMP_19_0.01$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.01, id = "result19", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.01, res_pa_GMP_19_0.01))

###########################################################  Prob = 0.05  #######################################################
prob <- 0.05
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.05", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.05 <- fit.classification(y=y, samples = samples, id = "result19_0.05", datapath = datapath, respath = respath, profile_name = profile_name,
                                           method = "DRW", pranking = "t-test", classifier = "rf",
                                           nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.05, file=file.path('data/model/res_pa_GMR_d_19_0.05.RData'))

summary(res_pa_GMR_d_19_0.05)
print(res_pa_GMR_d_19_0.05$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.05, id = "result19_0.05", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.05", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.05 <- fit.classification(y=y, samples = samples, id = "result19_0.05", datapath = datapath, respath = respath, profile_name = profile_name,
                                         method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.05, file=file.path('data/model/res_pa_GMP_19_0.05.RData'))

summary(res_pa_GMP_19_0.05)
print(res_pa_GMP_19_0.05$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.05, id = "result19_0.05", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.05, res_pa_GMP_19_0.05))

###########################################################  Prob = 0.1  #######################################################
prob <- 0.1
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.1", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.1 <- fit.classification(y=y, samples = samples, id = "result19_0.1", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.1, file=file.path('data/model/res_pa_GMR_d_19_0.1.RData'))

summary(res_pa_GMR_d_19_0.1)
print(res_pa_GMR_d_19_0.1$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.1, id = "result19_0.1", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.1", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.1 <- fit.classification(y=y, samples = samples, id = "result19_0.1", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.1, file=file.path('data/model/res_pa_GMP_19_0.1.RData'))

summary(res_pa_GMP_19_0.1)
print(res_pa_GMP_19_0.1$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.1, id = "result19_0.1", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.1, res_pa_GMP_19_0.1))

###########################################################  Prob = 0.2  #######################################################
prob <- 0.2
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.2", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.2 <- fit.classification(y=y, samples = samples, id = "result19_0.2", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.2, file=file.path('data/model/res_pa_GMR_d_19_0.2.RData'))

summary(res_pa_GMR_d_19_0.2)
print(res_pa_GMR_d_19_0.2$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.2, id = "result19_0.2", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.2", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.2 <- fit.classification(y=y, samples = samples, id = "result19_0.2", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.2, file=file.path('data/model/res_pa_GMP_19_0.2.RData'))

summary(res_pa_GMP_19_0.2)
print(res_pa_GMP_19_0.2$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.2, id = "result19_0.2", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.2, res_pa_GMP_19_0.2))

###########################################################  Prob = 0.3  #######################################################
prob <- 0.3
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.3", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.3 <- fit.classification(y=y, samples = samples, id = "result19_0.3", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.3, file=file.path('data/model/res_pa_GMR_d_19_0.3.RData'))

summary(res_pa_GMR_d_19_0.3)
print(res_pa_GMR_d_19_0.3$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.3, id = "result19_0.3", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.3", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.3 <- fit.classification(y=y, samples = samples, id = "result19_0.3", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.3, file=file.path('data/model/res_pa_GMP_19_0.3.RData'))

summary(res_pa_GMP_19_0.3)
print(res_pa_GMP_19_0.3$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.3, id = "result19_0.3", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.3, res_pa_GMP_19_0.3))

###########################################################  Prob = 0.4  #######################################################
prob <- 0.4
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.4", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.4 <- fit.classification(y=y, samples = samples, id = "result19_0.4", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.4, file=file.path('data/model/res_pa_GMR_d_19_0.4.RData'))

summary(res_pa_GMR_d_19_0.4)
print(res_pa_GMR_d_19_0.4$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.4, id = "result19_0.4", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.4", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.4 <- fit.classification(y=y, samples = samples, id = "result19_0.4", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.4, file=file.path('data/model/res_pa_GMP_19_0.4.RData'))

summary(res_pa_GMP_19_0.4)
print(res_pa_GMP_19_0.4$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.4, id = "result19_0.4", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.4, res_pa_GMP_19_0.4))

###########################################################  Prob = 0.5  #######################################################
prob <- 0.5
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.5", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.5 <- fit.classification(y=y, samples = samples, id = "result19_0.5", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.5, file=file.path('data/model/res_pa_GMR_d_19_0.5.RData'))

summary(res_pa_GMR_d_19_0.5)
print(res_pa_GMR_d_19_0.5$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.5, id = "result19_0.5", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.5", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.5 <- fit.classification(y=y, samples = samples, id = "result19_0.5", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.5, file=file.path('data/model/res_pa_GMP_19_0.5.RData'))

summary(res_pa_GMP_19_0.5)
print(res_pa_GMP_19_0.5$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.5, id = "result19_0.5", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.5, res_pa_GMP_19_0.5))

###########################################################  Prob = 0.6  #######################################################
prob <- 0.6
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.6", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.6 <- fit.classification(y=y, samples = samples, id = "result19_0.6", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.6, file=file.path('data/model/res_pa_GMR_d_19_0.6.RData'))

summary(res_pa_GMR_d_19_0.6)
print(res_pa_GMR_d_19_0.6$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.6, id = "result19_0.6", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.6", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.6 <- fit.classification(y=y, samples = samples, id = "result19_0.6", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.6, file=file.path('data/model/res_pa_GMP_19_0.6.RData'))

summary(res_pa_GMP_19_0.6)
print(res_pa_GMP_19_0.6$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.6, id = "result19_0.6", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.6, res_pa_GMP_19_0.6))

###########################################################  Prob = 0.7  #######################################################
prob <- 0.7
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.7", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.7 <- fit.classification(y=y, samples = samples, id = "result19_0.7", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.7, file=file.path('data/model/res_pa_GMR_d_19_0.7.RData'))

summary(res_pa_GMR_d_19_0.7)
print(res_pa_GMR_d_19_0.7$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.7, id = "result19_0.7", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.7", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.7 <- fit.classification(y=y, samples = samples, id = "result19_0.7", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.7, file=file.path('data/model/res_pa_GMP_19_0.7.RData'))

summary(res_pa_GMP_19_0.7)
print(res_pa_GMP_19_0.7$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.7, id = "result19_0.7", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.7, res_pa_GMP_19_0.7))


###########################################################  Prob = 0.8  #######################################################
prob <- 0.8
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.8", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.8 <- fit.classification(y=y, samples = samples, id = "result19_0.8", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.8, file=file.path('data/model/res_pa_GMR_d_19_0.8.RData'))

summary(res_pa_GMR_d_19_0.8)
print(res_pa_GMR_d_19_0.8$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.8, id = "result19_0.8", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.8", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.8 <- fit.classification(y=y, samples = samples, id = "result19_0.8", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.8, file=file.path('data/model/res_pa_GMP_19_0.8.RData'))

summary(res_pa_GMP_19_0.8)
print(res_pa_GMP_19_0.8$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.8, id = "result19_0.8", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.8, res_pa_GMP_19_0.8))

###########################################################  Prob = 0.9  #######################################################
prob <- 0.9
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.9", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.9 <- fit.classification(y=y, samples = samples, id = "result19_0.9", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.9, file=file.path('data/model/res_pa_GMR_d_19_0.9.RData'))

summary(res_pa_GMR_d_19_0.9)
print(res_pa_GMR_d_19_0.9$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.9, id = "result19_0.9", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.9", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.9 <- fit.classification(y=y, samples = samples, id = "result19_0.9", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.9, file=file.path('data/model/res_pa_GMP_19_0.9.RData'))

summary(res_pa_GMP_19_0.9)
print(res_pa_GMP_19_0.9$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.9, id = "result19_0.9", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.9, res_pa_GMP_19_0.9))


# plot
title <- c("Result 19")
xlabs <- c("GM", "GMR", "GMR_d_0.01", "GMP_0.01", "GMR_d_0.05", "GMP_0.05", "GMR_d_0.1", "GMP_0.1", "GMR_d_0.2", "GMP_0.2", "GMR_d_0.3", "GMP_0.3", "GMR_d_0.4", "GMP_0.4", "GMR_d_0.5", "GMP_0.5", "GMR_d_0.6", "GMP_0.6", "GMR_d_0.7", "GMP_0.7", "GMR_d_0.8", "GMP_0.8", "GMR_d_0.9", "GMP_0.9")
#res_models <- list(res_pa_GM_13, res_pa_GMR_13, res_pa_GMR_d_19_0.01, res_pa_GMP_19_0.01, res_pa_GMR_d_19_0.05, res_pa_GMP_19_0.05, res_pa_GMR_d_19_0.1, res_pa_GMP_19_0.1)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_RF_13$results$Accuracy[1])