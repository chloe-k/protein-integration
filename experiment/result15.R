# single gene profile on a feature data (updated in 2018/07/05)

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.

# DppiGraph_rdc.rda was used

# Classifier : rf(Random Forest)

#------------------------- mRNA expression gene profile -------------------------#
testStatistic <- c("t-test")
profile_name <- c('gf_rna(Entrez)')
x=list(rnaseq)

fit.iDRWPClass(x=x, y=y, testStatistic= testStatistic, profile_name = profile_name,
               datapath=datapath, respath = respath, method = "gf", samples = samples, id = "result15")

res_gf_G_15 <- fit.classification(y=y, samples = samples, id = "result15", datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "gf", pranking = "t-test", classifier = "rf",
                                  nFolds = 5, numTops=50, iter = 50)

save(res_gf_G_15, file=file.path('data/model/res_gf_G_15.RData'))

summary(res_gf_G_15)
print(res_gf_G_15$results)
print(res_gf_G_15$resample$Accuracy)


#------------------------- methylation expression gene profile -------------------------#
testStatistic <- c("t-test")
profile_name <- c('gf_meth(Entrez)')
x=list(imputed_methyl)

fit.iDRWPClass(x=x, y=y, testStatistic= testStatistic, profile_name = profile_name,
               datapath=datapath, respath = respath, method = "gf", samples = samples, id = "result15")

res_gf_M_15 <- fit.classification(y=y, samples = samples, id = "result15", datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "gf", pranking = "t-test", classifier = "rf",
                                  nFolds = 5, numTops=50, iter = 50)

save(res_gf_M_15, file=file.path('data/model/res_gf_M_15.RData'))

summary(res_gf_M_15)
print(res_gf_M_15$results)
print(res_gf_M_15$resample$Accuracy)


#------------------------- RPPA expression gene profile -------------------------#
testStatistic <- c("t-test")
profile_name <- c('gf_rppa(Entrez)')
x=list(rppa)

fit.iDRWPClass(x=x, y=y, testStatistic= testStatistic, profile_name = profile_name,
               datapath=datapath, respath = respath, method = "gf", samples = samples, id = "result15")

res_gf_P_15 <- fit.classification(y=y, samples = samples, id = "result15", datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "gf", pranking = "t-test", classifier = "rf",
                                  nFolds = 5, numTops=50, iter = 50)

save(res_gf_P_15, file=file.path('data/model/res_gf_P_15.RData'))

summary(res_gf_P_15)
print(res_gf_P_15$results)
print(res_gf_P_15$resample$Accuracy)


# plot
title <- c("Result 15")
xlabs <- c("G(gene)", "M(gene)", "P(gene)")
res_models <- list(res_gf_G_15, res_gf_M_15, res_gf_P_15)

perf_boxplot(title, xlabs, res_models, perf_min = 0.4, perf_max = 0.9, res_gf_G_15$results$Accuracy[1])