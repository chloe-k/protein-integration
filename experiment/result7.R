# single gene profile on a feature data (updated in 2018/06/05)

# Classifier : rf(Random Forest)

#------------------------- mRNA expression gene profile -------------------------#
testStatistic <- c("t-test")
profile_name <- c('gf_rna')
x=list(rnaseq)

fit.iDRWPClass(x=x, y=y, testStatistic= testStatistic, profile_name = profile_name,
               datapath=datapath, respath = respath, method = "gf", samples = samples)

res_gf_G_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                               method = "gf", pranking = "t-test", classifier = "rf",
                               nFolds = 5, numTops=50, iter = 50)

save(res_gf_G_RF, file=file.path('data/model/res_gf_G_RF.RData'))

summary(res_gf_G_RF)
print(res_gf_G_RF$results)
print(res_gf_G_RF$resample$Accuracy)


#------------------------- methylation expression gene profile -------------------------#
testStatistic <- c("t-test")
profile_name <- c('gf_meth')
x=list(imputed_methyl)

fit.iDRWPClass(x=x, y=y, testStatistic= testStatistic, profile_name = profile_name,
               datapath=datapath, respath = respath, method = "gf", samples = samples)

res_gf_M_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                               method = "gf", pranking = "t-test", classifier = "rf",
                               nFolds = 5, numTops=50, iter = 50)

save(res_gf_M_RF, file=file.path('data/model/res_gf_M_RF.RData'))

summary(res_gf_M_RF)
print(res_gf_M_RF$results)
print(res_gf_M_RF$resample$Accuracy)


#------------------------- RPPA expression gene profile -------------------------#
testStatistic <- c("t-test")
profile_name <- c('gf_rppa')
x=list(rppa)

fit.iDRWPClass(x=x, y=y, testStatistic= testStatistic, profile_name = profile_name,
               datapath=datapath, respath = respath, method = "gf", samples = samples)

res_gf_P_RF <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                               method = "gf", pranking = "t-test", classifier = "rf",
                               nFolds = 5, numTops=50, iter = 50)

save(res_gf_P_RF, file=file.path('data/model/res_gf_P_RF.RData'))

summary(res_gf_P_RF)
print(res_gf_P_RF$results)
print(res_gf_P_RF$resample$Accuracy)


# plot
xlabs <- c("G(gene)", "M(gene)", "P(gene)")
res_models <- list(res_gf_G_RF, res_gf_M_RF, res_gf_P_RF)

perf_boxplot(xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_gf_G_RF$results$Accuracy[1])