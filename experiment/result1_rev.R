# integrative DRW on combined feature data (updated in 2018/05/31)
# concat directed pathway graphs within each profile (GM & GMR & GMP)
# To check the effect of the edge direction, It has a reverse edge direction to result1's
# But, no difference between their result

# edge direction
# g -> m
# g -> p


#------------------------- RNAseq + Methyl(reverse edge direction) -------------------------#
 gm <- g %du% m
 testStatistic <- c("DESeq2", "t-test")
 profile_name <- c("rna(rev)", "meth(rev)")
 x=list(rnaseq, imputed_methyl)

 rev_res_pa_GM <- fit.iDRWPClass(x=x, y=y, globalGraph=gm,
                             testStatistic= testStatistic, profile_name = profile_name,
                             datapath = datapath, respath = respath, pathSet=pathSet,
                             method = "DRW", samples = samples, pranking = "t-test", mode = "GM",
                             nFolds = 5, iter = 50, AntiCorr=FALSE, DEBUG=TRUE)

 save(rev_res_pa_GM, file=file.path('data/model/rev_res_pa_GM.RData'))

 summary(rev_res_pa_GM)
 print(rev_res_pa_GM$results)
 print(rev_res_pa_GM$resample$Accuracy)

 write.SigFeatures(res_fit=rev_res_pa_GM, profile_name=profile_name, method="DRW", respath=respath)


 #------------------------- RNAseq + Methyl + RPPA(Pathway Graph)(reverse) -------------------------#
 gmr <- g %du% m %du% r
 testStatistic <- c("DESeq2", "t-test", "t-test")
 profile_name <- c("rna(rev)", "meth(rev)", "rppa(rev_Pathway_Graph)")
 x=list(rnaseq, imputed_methyl, rppa)

 rev_res_pa_GMR <- fit.iDRWPClass(x=x, y=y, globalGraph=gmr,
                              testStatistic= testStatistic, profile_name = profile_name,
                              datapath = datapath, respath = respath, pathSet=pathSet,
                              method = "DRW", samples = samples, pranking = "t-test", mode = "GMR",
                              nFolds = 5, iter = 50, AntiCorr=FALSE, DEBUG=TRUE)

 save(rev_res_pa_GMR, file=file.path('data/model/rev_res_pa_GMR.RData'))

 summary(rev_res_pa_GMR)
 print(rev_res_pa_GMR$results)
 print(rev_res_pa_GMR$resample$Accuracy)

 write.SigFeatures(res_fit=rev_res_pa_GMR, profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph)(reverse) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(rev)", "meth(rev)", "rppa(rev)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

rev_res_pa_GMP <- fit.iDRWPClass(x=x, y=y, globalGraph=gmp,
                              testStatistic= testStatistic, profile_name = profile_name,
                              datapath = datapath, respath = respath, pathSet=pathSet,
                              method = "DRW", samples = samples, pranking = "t-test", mode = "GMP",
                              nFolds = 5, iter = 50, AntiCorr=FALSE, DEBUG=TRUE)


save(rev_res_pa_GMP, file=file.path('data/model/rev_res_pa_GMP.RData'))

summary(rev_res_pa_GMP)
print(rev_res_pa_GMP$results)
print(rev_res_pa_GMP$resample$Accuracy)

write.SigFeatures(res_fit=rev_res_pa_GMP, profile_name=profile_name, method="DRW", respath=respath)

# plot
xlabs <- c("GM", "GMR", "GMP")
res_models <- list(rev_res_pa_GM, rev_res_pa_GMR, rev_res_pa_GMP)

perf_boxplot(xlabs, res_models, perf_min = 0.5, perf_max = 0.9, rev_res_pa_GM$results$Accuracy)