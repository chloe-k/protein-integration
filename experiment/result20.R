# integrative DRW on combined feature data (updated in 2018/07/06)
# concat directed pathway graphs within each profile (GM & GMR & GMR_d & GMP)

# In order to find optimized restart probability in DRW method.
# Grid search was performed in [0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.9, 0.8, 0.9]
# p=0.3(gamma) had used in before

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.

# Dppigraph(Entrez).rda was used

# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)


###########################################################  Gamma = 0.01  #######################################################
Gamma <- 0.01

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.01", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.01 <- fit.classification(y=y, samples = samples, id = "result20_0.01", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.01, file=file.path('data/model/res_pa_GM_20_0.01.RData'))

summary(res_pa_GM_20_0.01)
print(res_pa_GM_20_0.01$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.01, id = "result20_0.01", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.01", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.01 <- fit.classification(y=y, samples = samples, id = "result20_0.01", datapath = datapath, respath = respath, profile_name = profile_name,
                                         method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.01, file=file.path('data/model/res_pa_GMR_20_0.01.RData'))

summary(res_pa_GMR_20_0.01)
print(res_pa_GMR_20_0.01$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.01, id = "result20_0.01", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.01", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.01 <- fit.classification(y=y, samples = samples, id = "result20_0.01", datapath = datapath, respath = respath, profile_name = profile_name,
                                           method = "DRW", pranking = "t-test", classifier = "rf",
                                           nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.01, file=file.path('data/model/res_pa_GMR_d_20_0.01.RData'))

summary(res_pa_GMR_d_20_0.01)
print(res_pa_GMR_d_20_0.01$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.01, id = "result20_0.01", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.01", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.01 <- fit.classification(y=y, samples = samples, id = "result20_0.01", datapath = datapath, respath = respath, profile_name = profile_name,
                                         method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.01, file=file.path('data/model/res_pa_GMP_20_0.01.RData'))

summary(res_pa_GMP_20_0.01)
print(res_pa_GMP_20_0.01$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.01, id = "result20_0.01", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.01, res_pa_GMR_20_0.01, res_pa_GMR_d_20_0.01, res_pa_GMP_20_0.01))

###########################################################  Gamma = 0.05  #######################################################
Gamma <- 0.05

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.05", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.05 <- fit.classification(y=y, samples = samples, id = "result20_0.05", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.05, file=file.path('data/model/res_pa_GM_20_0.05.RData'))

summary(res_pa_GM_20_0.05)
print(res_pa_GM_20_0.05$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.05, id = "result20_0.05", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.05", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.05 <- fit.classification(y=y, samples = samples, id = "result20_0.05", datapath = datapath, respath = respath, profile_name = profile_name,
                                         method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.05, file=file.path('data/model/res_pa_GMR_20_0.05.RData'))

summary(res_pa_GMR_20_0.05)
print(res_pa_GMR_20_0.05$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.05, id = "result20_0.05", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.05", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.05 <- fit.classification(y=y, samples = samples, id = "result20_0.05", datapath = datapath, respath = respath, profile_name = profile_name,
                                           method = "DRW", pranking = "t-test", classifier = "rf",
                                           nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.05, file=file.path('data/model/res_pa_GMR_d_20_0.05.RData'))

summary(res_pa_GMR_d_20_0.05)
print(res_pa_GMR_d_20_0.05$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.05, id = "result20_0.05", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.05", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.05 <- fit.classification(y=y, samples = samples, id = "result20_0.05", datapath = datapath, respath = respath, profile_name = profile_name,
                                         method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.05, file=file.path('data/model/res_pa_GMP_20_0.05.RData'))

summary(res_pa_GMP_20_0.05)
print(res_pa_GMP_20_0.05$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.05, id = "result20_0.05", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.05, res_pa_GMR_20_0.05, res_pa_GMR_d_20_0.05, res_pa_GMP_20_0.05))

###########################################################  Gamma = 0.1  #######################################################
Gamma <- 0.1

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.1", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.1 <- fit.classification(y=y, samples = samples, id = "result20_0.1", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.1, file=file.path('data/model/res_pa_GM_20_0.1.RData'))

summary(res_pa_GM_20_0.1)
print(res_pa_GM_20_0.1$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.1, id = "result20_0.1", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.1", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.1 <- fit.classification(y=y, samples = samples, id = "result20_0.1", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.1, file=file.path('data/model/res_pa_GMR_20_0.1.RData'))

summary(res_pa_GMR_20_0.1)
print(res_pa_GMR_20_0.1$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.1, id = "result20_0.1", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.1", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.1 <- fit.classification(y=y, samples = samples, id = "result20_0.1", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.1, file=file.path('data/model/res_pa_GMR_d_20_0.1.RData'))

summary(res_pa_GMR_d_20_0.1)
print(res_pa_GMR_d_20_0.1$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.1, id = "result20_0.1", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.1", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.1 <- fit.classification(y=y, samples = samples, id = "result20_0.1", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.1, file=file.path('data/model/res_pa_GMP_20_0.1.RData'))

summary(res_pa_GMP_20_0.1)
print(res_pa_GMP_20_0.1$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.1, id = "result20_0.1", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.1, res_pa_GMR_20_0.1, res_pa_GMR_d_20_0.1, res_pa_GMP_20_0.1))

###########################################################  Gamma = 0.2  #######################################################
Gamma <- 0.2

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.2", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.2 <- fit.classification(y=y, samples = samples, id = "result20_0.2", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.2, file=file.path('data/model/res_pa_GM_20_0.2.RData'))

summary(res_pa_GM_20_0.2)
print(res_pa_GM_20_0.2$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.2, id = "result20_0.2", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.2", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.2 <- fit.classification(y=y, samples = samples, id = "result20_0.2", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.2, file=file.path('data/model/res_pa_GMR_20_0.2.RData'))

summary(res_pa_GMR_20_0.2)
print(res_pa_GMR_20_0.2$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.2, id = "result20_0.2", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.2", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.2 <- fit.classification(y=y, samples = samples, id = "result20_0.2", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.2, file=file.path('data/model/res_pa_GMR_d_20_0.2.RData'))

summary(res_pa_GMR_d_20_0.2)
print(res_pa_GMR_d_20_0.2$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.2, id = "result20_0.2", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.2", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.2 <- fit.classification(y=y, samples = samples, id = "result20_0.2", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.2, file=file.path('data/model/res_pa_GMP_20_0.2.RData'))

summary(res_pa_GMP_20_0.2)
print(res_pa_GMP_20_0.2$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.2, id = "result20_0.2", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.2, res_pa_GMR_20_0.2, res_pa_GMR_d_20_0.2, res_pa_GMP_20_0.2))

###########################################################  Gamma = 0.3  #######################################################
Gamma <- 0.3

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.3", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.3 <- fit.classification(y=y, samples = samples, id = "result20_0.3", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.3, file=file.path('data/model/res_pa_GM_20_0.3.RData'))

summary(res_pa_GM_20_0.3)
print(res_pa_GM_20_0.3$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.3, id = "result20_0.3", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.3", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.3 <- fit.classification(y=y, samples = samples, id = "result20_0.3", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.3, file=file.path('data/model/res_pa_GMR_20_0.3.RData'))

summary(res_pa_GMR_20_0.3)
print(res_pa_GMR_20_0.3$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.3, id = "result20_0.3", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.3", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.3 <- fit.classification(y=y, samples = samples, id = "result20_0.3", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.3, file=file.path('data/model/res_pa_GMR_d_20_0.3.RData'))

summary(res_pa_GMR_d_20_0.3)
print(res_pa_GMR_d_20_0.3$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.3, id = "result20_0.3", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.3", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.3 <- fit.classification(y=y, samples = samples, id = "result20_0.3", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.3, file=file.path('data/model/res_pa_GMP_20_0.3.RData'))

summary(res_pa_GMP_20_0.3)
print(res_pa_GMP_20_0.3$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.3, id = "result20_0.3", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.3, res_pa_GMR_20_0.3, res_pa_GMR_d_20_0.3, res_pa_GMP_20_0.3))


###########################################################  Gamma = 0.4  #######################################################
Gamma <- 0.4

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.4", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.4 <- fit.classification(y=y, samples = samples, id = "result20_0.4", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.4, file=file.path('data/model/res_pa_GM_20_0.4.RData'))

summary(res_pa_GM_20_0.4)
print(res_pa_GM_20_0.4$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.4, id = "result20_0.4", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.4", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.4 <- fit.classification(y=y, samples = samples, id = "result20_0.4", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.4, file=file.path('data/model/res_pa_GMR_20_0.4.RData'))

summary(res_pa_GMR_20_0.4)
print(res_pa_GMR_20_0.4$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.4, id = "result20_0.4", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.4", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.4 <- fit.classification(y=y, samples = samples, id = "result20_0.4", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.4, file=file.path('data/model/res_pa_GMR_d_20_0.4.RData'))

summary(res_pa_GMR_d_20_0.4)
print(res_pa_GMR_d_20_0.4$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.4, id = "result20_0.4", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.4", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.4 <- fit.classification(y=y, samples = samples, id = "result20_0.4", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.4, file=file.path('data/model/res_pa_GMP_20_0.4.RData'))

summary(res_pa_GMP_20_0.4)
print(res_pa_GMP_20_0.4$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.4, id = "result20_0.4", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.4, res_pa_GMR_20_0.4, res_pa_GMR_d_20_0.4, res_pa_GMP_20_0.4))

###########################################################  Gamma = 0.5  #######################################################
Gamma <- 0.5

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.5", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.5 <- fit.classification(y=y, samples = samples, id = "result20_0.5", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.5, file=file.path('data/model/res_pa_GM_20_0.5.RData'))

summary(res_pa_GM_20_0.5)
print(res_pa_GM_20_0.5$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.5, id = "result20_0.5", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.5", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.5 <- fit.classification(y=y, samples = samples, id = "result20_0.5", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.5, file=file.path('data/model/res_pa_GMR_20_0.5.RData'))

summary(res_pa_GMR_20_0.5)
print(res_pa_GMR_20_0.5$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.5, id = "result20_0.5", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.5", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.5 <- fit.classification(y=y, samples = samples, id = "result20_0.5", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.5, file=file.path('data/model/res_pa_GMR_d_20_0.5.RData'))

summary(res_pa_GMR_d_20_0.5)
print(res_pa_GMR_d_20_0.5$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.5, id = "result20_0.5", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.5", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.5 <- fit.classification(y=y, samples = samples, id = "result20_0.5", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.5, file=file.path('data/model/res_pa_GMP_20_0.5.RData'))

summary(res_pa_GMP_20_0.5)
print(res_pa_GMP_20_0.5$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.5, id = "result20_0.5", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.5, res_pa_GMR_20_0.5, res_pa_GMR_d_20_0.5, res_pa_GMP_20_0.5))


###########################################################  Gamma = 0.6  #######################################################
Gamma <- 0.6

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.6", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.6 <- fit.classification(y=y, samples = samples, id = "result20_0.6", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.6, file=file.path('data/model/res_pa_GM_20_0.6.RData'))

summary(res_pa_GM_20_0.6)
print(res_pa_GM_20_0.6$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.6, id = "result20_0.6", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.6", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.6 <- fit.classification(y=y, samples = samples, id = "result20_0.6", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.6, file=file.path('data/model/res_pa_GMR_20_0.6.RData'))

summary(res_pa_GMR_20_0.6)
print(res_pa_GMR_20_0.6$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.6, id = "result20_0.6", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.6", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.6 <- fit.classification(y=y, samples = samples, id = "result20_0.6", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.6, file=file.path('data/model/res_pa_GMR_d_20_0.6.RData'))

summary(res_pa_GMR_d_20_0.6)
print(res_pa_GMR_d_20_0.6$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.6, id = "result20_0.6", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.6", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.6 <- fit.classification(y=y, samples = samples, id = "result20_0.6", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.6, file=file.path('data/model/res_pa_GMP_20_0.6.RData'))

summary(res_pa_GMP_20_0.6)
print(res_pa_GMP_20_0.6$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.6, id = "result20_0.6", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.6, res_pa_GMR_20_0.6, res_pa_GMR_d_20_0.6, res_pa_GMP_20_0.6))

###########################################################  Gamma = 0.7  #######################################################
Gamma <- 0.7

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.7", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.7 <- fit.classification(y=y, samples = samples, id = "result20_0.7", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.7, file=file.path('data/model/res_pa_GM_20_0.7.RData'))

summary(res_pa_GM_20_0.7)
print(res_pa_GM_20_0.7$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.7, id = "result20_0.7", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.7", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.7 <- fit.classification(y=y, samples = samples, id = "result20_0.7", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.7, file=file.path('data/model/res_pa_GMR_20_0.7.RData'))

summary(res_pa_GMR_20_0.7)
print(res_pa_GMR_20_0.7$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.7, id = "result20_0.7", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.7", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.7 <- fit.classification(y=y, samples = samples, id = "result20_0.7", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.7, file=file.path('data/model/res_pa_GMR_d_20_0.7.RData'))

summary(res_pa_GMR_d_20_0.7)
print(res_pa_GMR_d_20_0.7$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.7, id = "result20_0.7", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.7", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.7 <- fit.classification(y=y, samples = samples, id = "result20_0.7", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.7, file=file.path('data/model/res_pa_GMP_20_0.7.RData'))

summary(res_pa_GMP_20_0.7)
print(res_pa_GMP_20_0.7$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.7, id = "result20_0.7", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.7, res_pa_GMR_20_0.7, res_pa_GMR_d_20_0.7, res_pa_GMP_20_0.7))

###########################################################  Gamma = 0.8  #######################################################
Gamma <- 0.8

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.8", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.8 <- fit.classification(y=y, samples = samples, id = "result20_0.8", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.8, file=file.path('data/model/res_pa_GM_20_0.8.RData'))

summary(res_pa_GM_20_0.8)
print(res_pa_GM_20_0.8$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.8, id = "result20_0.8", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.8", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.8 <- fit.classification(y=y, samples = samples, id = "result20_0.8", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.8, file=file.path('data/model/res_pa_GMR_20_0.8.RData'))

summary(res_pa_GMR_20_0.8)
print(res_pa_GMR_20_0.8$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.8, id = "result20_0.8", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.8", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.8 <- fit.classification(y=y, samples = samples, id = "result20_0.8", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.8, file=file.path('data/model/res_pa_GMR_d_20_0.8.RData'))

summary(res_pa_GMR_d_20_0.8)
print(res_pa_GMR_d_20_0.8$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.8, id = "result20_0.8", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.8", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.8 <- fit.classification(y=y, samples = samples, id = "result20_0.8", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.8, file=file.path('data/model/res_pa_GMP_20_0.8.RData'))

summary(res_pa_GMP_20_0.8)
print(res_pa_GMP_20_0.8$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.8, id = "result20_0.8", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.8, res_pa_GMR_20_0.8, res_pa_GMR_d_20_0.8, res_pa_GMP_20_0.8))


###########################################################  Gamma = 0.9  #######################################################
Gamma <- 0.9

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.9", Gamma = Gamma,
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20_0.9 <- fit.classification(y=y, samples = samples, id = "result20_0.9", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_20_0.9, file=file.path('data/model/res_pa_GM_20_0.9.RData'))

summary(res_pa_GM_20_0.9)
print(res_pa_GM_20_0.9$results)

write.SigFeatures(res_fit=res_pa_GM_20_0.9, id = "result20_0.9", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.9", Gamma = Gamma,
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20_0.9 <- fit.classification(y=y, samples = samples, id = "result20_0.9", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_20_0.9, file=file.path('data/model/res_pa_GMR_20_0.9.RData'))

summary(res_pa_GMR_20_0.9)
print(res_pa_GMR_20_0.9$results)

write.SigFeatures(res_fit=res_pa_GMR_20_0.9, id = "result20_0.9", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.9", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20_0.9 <- fit.classification(y=y, samples = samples, id = "result20_0.9", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_20_0.9, file=file.path('data/model/res_pa_GMR_d_20_0.9.RData'))

summary(res_pa_GMR_d_20_0.9)
print(res_pa_GMR_d_20_0.9$results)

write.SigFeatures(res_fit=res_pa_GMR_d_20_0.9, id = "result20_0.9", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result20_0.9", Gamma = Gamma,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20_0.9 <- fit.classification(y=y, samples = samples, id = "result20_0.9", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_20_0.9, file=file.path('data/model/res_pa_GMP_20_0.9.RData'))

summary(res_pa_GMP_20_0.9)
print(res_pa_GMP_20_0.9$results)

write.SigFeatures(res_fit=res_pa_GMP_20_0.9, id = "result20_0.9", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20_0.9, res_pa_GMR_20_0.9, res_pa_GMR_d_20_0.9, res_pa_GMP_20_0.9))


# plot
title <- c("Result 20")
xlabs <- c("GM_0.01", "GMR_0.01", "GMR_d_0.01", "GMP_0.01", "GM_0.05", "GMR_0.05", "GMR_d_0.05", "GMP_0.05", "GM_0.1", "GMR_0.1", "GMR_d_0.1", "GMP_0.1", "GM_0.2", "GMR_0.2", "GMR_d_0.2", "GMP_0.2", "GM_0.3", "GMR_0.3", "GMR_d_0.3", "GMP_0.3","GM_0.4", "GMR_0.4", "GMR_d_0.4", "GMP_0.4", "GM_0.5", "GMR_0.5", "GMR_d_0.5", "GMP_0.5", "GM_0.6", "GMR_0.6", "GMR_d_0.6", "GMP_0.6", "GM_0.7", "GMR_0.7", "GMR_d_0.7", "GMP_0.7", "GM_0.8", "GMR_0.8", "GMR_d_0.8", "GMP_0.8", "GM_0.9", "GMR_0.9", "GMR_d_0.9", "GMP_0.9")

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, mean(res_pa_GM_20_0.01$results$Accuracy))

# plot for each models
res_gm <- list()
res_gmr <- list()
res_gmr_d <- list()
res_gmp <- list()
for(i in 1:length(res_models)){
  if(i %% 4 == 1){
    res_gm <- c(res_gm, list(res_models[[i]]))
  }
  if(i %% 4 == 2){
    res_gmr <- c(res_gmr, list(res_models[[i]]))
  }
  if(i %% 4 == 3){
    res_gmr_d <- c(res_gmr_d, list(res_models[[i]]))
  }
  if(i %% 4 == 0){
    res_gmp <- c(res_gmp, list(res_models[[i]]))
  }
}

# plot
title <- c("Result GM")
xlabs <- c("GM_0.01", "GM_0.05", "GM_0.1", "GM_0.2", "GM_0.3", "GM_0.4", "GM_0.5", "GM_0.6", "GM_0.7", "GM_0.8", "GM_0.9")
perf_boxplot(title, xlabs, res_gm, perf_min = 0.5, perf_max = 0.9, mean(sapply(X = res_gm, FUN = function(x){x$results$Accuracy})))

# plot
title <- c("Result GMR")
xlabs <- c("GMR_0.01", "GMR_0.05", "GMR_0.1", "GMR_0.2", "GMR_0.3", "GMR_0.4", "GMR_0.5", "GMR_0.6", "GMR_0.7", "GMR_0.8", "GMR_0.9")
perf_boxplot(title, xlabs, res_gmr, perf_min = 0.5, perf_max = 0.9, mean(sapply(X = res_gmr, FUN = function(x){x$results$Accuracy})))

# plot
title <- c("Result GMR_d")
xlabs <- c("GMR_d_0.01", "GMR_d_0.05", "GMR_d_0.1", "GMR_d_0.2", "GMR_d_0.3", "GMR_d_0.4", "GMR_d_0.5", "GMR_d_0.6", "GMR_d_0.7", "GMR_d_0.8", "GMR_d_0.9")
perf_boxplot(title, xlabs, res_gmr_d, perf_min = 0.5, perf_max = 0.9, mean(sapply(X = res_gmr_d, FUN = function(x){x$results$Accuracy})))

# plot
title <- c("Result GMP")
xlabs <- c("GMP_0.01", "GMP_0.05", "GMP_0.1", "GMP_0.2", "GMP_0.3", "GMP_0.4", "GMP_0.5", "GMP_0.6", "GMP_0.7", "GMP_0.8", "GMP_0.9")
perf_boxplot(title, xlabs, res_gmp, perf_min = 0.5, perf_max = 0.9, mean(sapply(X = res_gmp, FUN = function(x){x$results$Accuracy})))