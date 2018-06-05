# single gene profile on a feature data (updated in 2018/05/31)

# Classifier : glm(Generalized Linear Model)

#------------------------- mRNA expression gene profile -------------------------#
fit.iDRWPClass(x=list(rnaseq), y=y, testStatistic= c('t-test'), profile_name = c('gf_rna'),
                             datapath=datapath, respath = respath, method = "gf", samples = samples)

res_gf_G <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                               method = "DRW", pranking = "t-test", classifier = "glm",
                               nFolds = 5, numTops=50, iter = 50)

save(res_gf_G, file=file.path('data/model/res_gf_G.RData'))

summary(res_gf_G)
print(res_gf_G$results)
print(res_gf_G$resample$Accuracy)


#------------------------- methylation expression gene profile -------------------------#
fit.iDRWPClass(x=list(imputed_methyl), y=y, testStatistic= c('t-test'), profile_name = c('gf_meth'),
                           datapath=datapath, respath = respath, method = "gf", samples = samples)

res_gf_M <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                               method = "DRW", pranking = "t-test", classifier = "glm",
                               nFolds = 5, numTops=50, iter = 50)

save(res_gf_M, file=file.path('data/model/res_gf_M.RData'))

summary(res_gf_M)
print(res_gf_M$results)
print(res_gf_M$resample$Accuracy)


#------------------------- RPPA expression gene profile -------------------------#
fit.iDRWPClass(x=list(rppa), y=y, testStatistic= c('t-test'), profile_name = c('gf_rppa'),
                           datapath=datapath, respath = respath, method = "gf", samples = samples)

res_gf_P <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                               method = "DRW", pranking = "t-test", classifier = "glm",
                               nFolds = 5, numTops=50, iter = 50)

save(res_gf_P, file=file.path('data/model/res_gf_P.RData'))

summary(res_gf_P)
print(res_gf_P$results)
print(res_gf_P$resample$Accuracy)


# plot
xlabs <- c("G(gene)", "M(gene)", "P(gene)")
res_models <- list(res_gf_G, res_gf_M, res_gf_P)

perf_boxplot(xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_gf_G$results$Accuracy)