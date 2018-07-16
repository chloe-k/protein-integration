# integrative DRW on combined feature data (updated in 2018/07/16)
# concat directed pathway graphs within each profile (G & M & R & GM & GR & MP & GMR)

# p=0.4, g=0.2 was used for restart probability
# When pathway activity score was calculated, each type of weight was used.
# This experiment is based on GMR_d model
# G -> perform GMR_d model and only 'g' type weight is included for pathwway activity score calculation 

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.


# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)

################################## Result 23 ############################################################

res_models <- list()


################################################### Result23 #################################################


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) - G -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result23_G", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_G_23 <- fit.classification(y=y, samples = samples, id = "result23_G", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      numTops=50)

save(res_pa_G_23, file=file.path('data/model/res_pa_G_23.RData'))

summary(res_pa_G_23)

write.SigFeatures(res_fit=res_pa_G_23, id = "result23_G", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_G_23))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) - M -------------------------#
sapply(file.path("utils",list.files("utils", pattern="*.R")),source)
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result23_M", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_M_23 <- fit.classification(y=y, samples = samples, id = "result23_M", datapath = datapath, respath = respath,
                                  profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                  numTops=50)

save(res_pa_M_23, file=file.path('data/model/res_pa_M_23.RData'))

summary(res_pa_M_23)

write.SigFeatures(res_fit=res_pa_M_23, id = "result23_M", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_M_23))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) - R -------------------------#
sapply(file.path("utils",list.files("utils", pattern="*.R")),source)
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result23_R", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_R_23 <- fit.classification(y=y, samples = samples, id = "result23_R", datapath = datapath, respath = respath,
                                  profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                  numTops=50)

save(res_pa_R_23, file=file.path('data/model/res_pa_R_23.RData'))

summary(res_pa_R_23)

write.SigFeatures(res_fit=res_pa_R_23, id = "result23_R", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_R_23))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) - GM -------------------------#
sapply(file.path("utils",list.files("utils", pattern="*.R")),source)
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result23_GM", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_23 <- fit.classification(y=y, samples = samples, id = "result23_GM", datapath = datapath, respath = respath,
                                  profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                  numTops=50)

save(res_pa_GM_23, file=file.path('data/model/res_pa_GM_23.RData'))

summary(res_pa_GM_23)

write.SigFeatures(res_fit=res_pa_GM_23, id = "result23_GM", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_23))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) - GR -------------------------#
sapply(file.path("utils",list.files("utils", pattern="*.R")),source)
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result23_GR", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GR_23 <- fit.classification(y=y, samples = samples, id = "result23_GR", datapath = datapath, respath = respath,
                                   profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                   numTops=50)

save(res_pa_GR_23, file=file.path('data/model/res_pa_GR_23.RData'))

summary(res_pa_GR_23)

write.SigFeatures(res_fit=res_pa_GR_23, id = "result23_GR", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GR_23))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) - MR -------------------------#
sapply(file.path("utils",list.files("utils", pattern="*.R")),source)
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result23_MR", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_MR_23 <- fit.classification(y=y, samples = samples, id = "result23_MR", datapath = datapath, respath = respath,
                                   profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                   numTops=50)

save(res_pa_MR_23, file=file.path('data/model/res_pa_MR_23.RData'))

summary(res_pa_MR_23)

write.SigFeatures(res_fit=res_pa_MR_23, id = "result23_MR", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_MR_23))



#########################################################################################################################################
# plot

title <- c("Result 23")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_models, FUN = function(x){mean(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){mean(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.02, perf_max = perf_max+0.02)
