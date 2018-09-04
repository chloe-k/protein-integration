# single gene profile on a feature data (updated in 2018/06/25)

# All gene symbols are converted to Entrez gene id
# LOOCV was performed.

# Dppigraph(Entrez).rda was used

# Classifier : rf(Random Forest)

testStatistic <- c("t-test")
profile_name <- c('gf_rppa(Entrez)')
x=list(rppa)

fit.iDRWPClass(x=x, y=y, testStatistic= testStatistic, profile_name = profile_name, mode = NULL, id = "result12_P",
               datapath=datapath, respath = respath, method = "gf", samples = samples)

res_gf_P_12 <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "gf", pranking = "t-test", classifier = "rf", id = "result12_P",
                                  nFolds = 5, numTops=50, iter = 10)

save(res_gf_P_12, file=file.path('data/model/res_gf_P_12.RData'))

#------------------------------------
testStatistic <- c("t-test")
profile_name <- c('gf_meth(Entrez)')
x=list(imputed_methyl)

fit.iDRWPClass(x=x, y=y, testStatistic= testStatistic, profile_name = profile_name, mode = NULL, id = "result12_M",
               datapath=datapath, respath = respath, method = "gf", samples = samples)

res_gf_M_12 <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "gf", pranking = "t-test", classifier = "rf", id = "result12_M",
                                  nFolds = 5, numTops=50, iter = 10)

save(res_gf_M_12, file=file.path('data/model/res_gf_M_12.RData'))

#----------------------
testStatistic <- c('DESeq2')
profile_name <- c('gf_rna(Entrez)')
x=list(rnaseq)

fit.iDRWPClass(x=x, y=y, testStatistic= testStatistic, profile_name = profile_name, mode = NULL, id = "result12_G",
               datapath=datapath, respath = respath, method = "gf", samples = samples)

res_gf_G_12 <- fit.classification(y=y, samples = samples, datapath = datapath, respath = respath, profile_name = profile_name,
                                  method = "gf", pranking = "t-test", classifier = "rf", id = "result12_G",
                                  nFolds = 5, numTops=50, iter = 10)

save(res_gf_G_12, file=file.path('data/model/res_gf_G_12.RData'))


# plot
title <- c("Result 12")
xlabs <- c("G(gene)", "M(gene)", "P(protein)")
res_models <- list(res_gf_G_12, res_gf_M_12, res_gf_P_12)

perf_boxplot(title, xlabs, res_models, perf_min = 0.4, perf_max = 0.9, res_gf_G_12$results$Accuracy[1])