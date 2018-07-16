# integrative DRW on combined feature data (updated in 2018/07/16)
# concat directed pathway graphs within each profile (GM & GMR & GMR_d & GMP & GMR_d_str & GMP_str)

# p=0.4, g=0.2 was used for restart probability
# Dppigraph(Entrez).rda was used in GM, GMR, GMR_d, GMP 
# DppiGraph_W_str.rda was used in GMR_d_str, GMP_str

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.


# edge direction
# m -> g
# p -> g

# Classifier : rf(Random Forest)

################################## Result 20 ############################################################

res_models <- list()


################################################### Result20: prob = 0.4, Gamma = 0.2  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result20_GM", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_20 <- fit.classification(y=y, samples = samples, id = "result20_GM", datapath = datapath, respath = respath,
                                   profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                   numTops=50)

save(res_pa_GM_20, file=file.path('data/model/res_pa_GM_20.RData'))

summary(res_pa_GM_20)

write.SigFeatures(res_fit=res_pa_GM_20, id = "result20_GM", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_20))

#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result20_GMR", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_20 <- fit.classification(y=y, samples = samples, id = "result20_GMR", datapath = datapath, respath = respath,
                                    profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                    numTops=50)

save(res_pa_GMR_20, file=file.path('data/model/res_pa_GMR_20.RData'))

summary(res_pa_GMR_20)

write.SigFeatures(res_fit=res_pa_GMR_20, id = "result20_GMR", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_20))

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result20_GMR_d", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_20 <- fit.classification(y=y, samples = samples, id = "result20_GMR_d", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      numTops=50)

save(res_pa_GMR_d_20, file=file.path('data/model/res_pa_GMR_d_20.RData'))

summary(res_pa_GMR_d_20)

write.SigFeatures(res_fit=res_pa_GMR_d_20, id = "result20_GMR_d", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_20))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result20_GMP", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_20 <- fit.classification(y=y, samples = samples, id = "result20_GMP", datapath = datapath, respath = respath,
                                    profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                    numTops=50)


save(res_pa_GMP_20, file=file.path('data/model/res_pa_GMP_20.RData'))

summary(res_pa_GMP_20)

write.SigFeatures(res_fit=res_pa_GMP_20, id = "result20_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMP_20))


#------------------------- RNAseq + Methyl + RPPA(diffused weighted STRING Pathway Graph) -------------------------#
dppipath <- file.path(datapath, 'directGraph_W_str.rda')
load(file.path(dppipath))

r <- directGraph
V(r)$name <-paste("p",V(r)$name,sep="")

gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result20_GMR_d_w", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_w_20 <- fit.classification(y=y, samples = samples, id = "result20_GMR_d_w", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      numTops=50)

save(res_pa_GMR_d_w_20, file=file.path('data/model/res_pa_GMR_d_w_20.RData'))

summary(res_pa_GMR_d_w_20)

write.SigFeatures(res_fit=res_pa_GMR_d_w_20, id = "result20_GMR_d_w", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_w_20))


#------------------------- RNAseq + Methyl + RPPA(weighted STRING PPI Graph) -------------------------#
dppipath <- file.path(datapath, 'DppiGraph_W_str.rda')
load(file.path(dppipath))

p <- DppiGraph
V(p)$name <-paste("p",V(p)$name,sep="")

testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result20_GMP_w", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_w_20 <- fit.classification(y=y, samples = samples, id = "result20_GMP_w", datapath = datapath, respath = respath,
                                    profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                    numTops=50)


save(res_pa_GMP_w_20, file=file.path('data/model/res_pa_GMP_w_20.RData'))

summary(res_pa_GMP_w_20)

write.SigFeatures(res_fit=res_pa_GMP_w_20, id = "result20_GMP_w", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMP_w_20))


#########################################################################################################################################
# plot

title <- c("Result 20")
xlabs <- c("GM", "GMR", "GMR_d", "GMP", "GMR_d_str", "GMP_str")

perf_min <- min(sapply(X = res_models, FUN = function(x){mean(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){mean(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.02, perf_max = perf_max+0.02)
