
########################### package load ############################
# source("https://bioconductor.org/biocLite.R")

library(KEGGgraph)
library(igraph)
library(ggplot2)
library(annotate)
library(org.Hs.eg.db)
library(diffusr)
library(DESeq2)
library(Matrix)
library(stringr)
library(caret)
library(e1071)
library(randomForest)
library(KEGG.db)
library(KEGGREST)
library(biomaRt)

######################################################################

sapply(file.path("utils",list.files("utils", pattern="*.R")),source)

# make directory
datapath <- file.path('data')
if(!dir.exists(datapath)) dir.create(datapath)

gdacpath <- file.path(datapath, 'BRCA_GDAC')

respath <- file.path('result')
if(!dir.exists(respath)) dir.create(respath)


# read RNAseq, Methylation data, RPPA data
#data_all_path <- file.path(datapath, "data.RData")
data_all_path <- file.path(datapath, "Entrez_data.RData")
if(!file.exists(data_all_path)) {
  year <- 3
  read_data(year, datapath)
}
load(data_all_path)


# read rda data
# graphpath <- file.path(datapath,'directGraph.rda')
# pathSetpath <- file.path(datapath,'pathSet.rda')
graphpath <- file.path(datapath,'directGraph(Entrez).rda')
pathSetpath <- file.path(datapath,'pathSet(Entrez).rda')
if(!(file.exists(graphpath) && file.exists(pathSetpath))) {
  print('directGraph and pathSet do not exist, now creating directGraph and pathSet is start')
  cons_KEGGgraph(datapath)
}
load(file.path(graphpath))
load(file.path(pathSetpath))



dppipath <- file.path(datapath, 'DppiGraph(Entrez).rda')  # directed edge PPI
# dppipath <- file.path(datapath, 'DppiGraph_rdc.rda')  # directed edge PPI
# dppipath <- file.path(datapath, 'DppiGraph_W_str.rda')

if(!file.exists(dppipath)) {
  print('ppiGraph does not exist, now creating ppiGraph start')
  cons_ppi(datapath, gdacpath, rppa)
}
load(file.path(dppipath))


# directed pathway graph provided in DRWPClass
g <- directGraph 
V(g)$name <- paste("g",V(g)$name,sep="")

m <- directGraph
V(m)$name <-paste("m",V(m)$name,sep="")

r <- directGraph
V(r)$name <-paste("p",V(r)$name,sep="")

p <- DppiGraph
V(p)$name <-paste("p",V(p)$name,sep="")

y=list(good_samples, poor_samples)
#---------------DRW---------------#
# RNAseq pathway profile


# Methylation pathway profile


# RPPA pathway profile


#----------------------------------------iDRW-----------------------------------------------------------#

################################## Result 18_all ############################################################

res_gm <- list()

################################################### Result18_1: prob = 0.001, Gamma = 0  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_1_GM", prob = 0.001, Gamma = 0, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_1 <- fit.classification(y=y, samples = samples, id = "result18_1_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_1, file=file.path('data/model/res_pa_GM_18_1.RData'))

res_gm <- c(res_gm, list(res_pa_GM_18_1))

################################################### Result18_2: prob = 0.001, Gamma = 0.2  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_2_GM", prob = 0.001, Gamma = 0.2, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_2 <- fit.classification(y=y, samples = samples, id = "result18_2_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_2, file=file.path('data/model/res_pa_GM_18_2.RData'))

res_gm <- c(res_gm, list(res_pa_GM_18_2))

################################################### Result18_3: prob = 0.001, Gamma = 0.4  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_3_GM", prob = 0.001, Gamma = 0.4, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_3 <- fit.classification(y=y, samples = samples, id = "result18_3_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_3, file=file.path('data/model/res_pa_GM_18_3.RData'))

res_gm <- c(res_gm, list(res_pa_GM_18_3))

################################################### Result18_4: prob = 0.001, Gamma = 0.6  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_4_GM", prob = 0.001, Gamma = 0.6, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_4 <- fit.classification(y=y, samples = samples, id = "result18_4_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_4, file=file.path('data/model/res_pa_GM_18_4.RData'))

res_gm <- c(res_gm, list(res_pa_GM_18_4))

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_4.5_GM", prob = 0.001, Gamma = 0.7, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_4.5 <- fit.classification(y=y, samples = samples, id = "result18_4.5_GM", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_4.5, file=file.path('data/model/res_pa_GM_18_4.5.RData'))

res_gm <- c(res_gm, list(res_pa_GM_18_4.5))

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_0.75_GM", prob = 0.001, Gamma = 0.75, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_0.75 <- fit.classification(y=y, samples = samples, id = "result18_0.75_GM", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_0.75, file=file.path('data/model/res_pa_GM_18_0.75.RData'))

res_gm <- c(res_gm, list(res_pa_GM_18_0.75))

################################################### Result18_5: prob = 0.001, Gamma = 0.8  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_5_GM", prob = 0.001, Gamma = 0.8, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_5 <- fit.classification(y=y, samples = samples, id = "result18_5_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_5, file=file.path('data/model/res_pa_GM_18_5.RData'))

res_gm <- c(res_gm, list(res_pa_GM_18_5))


#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_0.85_GM", prob = 0.001, Gamma = 0.85, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_0.85 <- fit.classification(y=y, samples = samples, id = "result18_0.85_GM", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_0.85, file=file.path('data/model/res_pa_GM_18_0.85.RData'))

res_gm <- c(res_gm, list(res_pa_GM_18_0.85))

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_0.9_GM", prob = 0.001, Gamma = 0.9, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_0.9 <- fit.classification(y=y, samples = samples, id = "result18_0.9_GM", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_0.9, file=file.path('data/model/res_pa_GM_18_0.9.RData'))

res_gm <- c(res_gm, list(res_pa_GM_18_0.9))


#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_0.95_GM", prob = 0.001, Gamma = 0.95, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_0.95 <- fit.classification(y=y, samples = samples, id = "result18_0.95_GM", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_0.95, file=file.path('data/model/res_pa_GM_18_0.95.RData'))

res_gm <- c(res_gm, list(res_pa_GM_18_0.95))


#########################################################################################################################################
# plot

# GM
title <- c("Result 18_GM")
xlabs <- c("[g=0]", "[g=0.2]", "[g=0.4]", "[g=0.6]", "[g=0.7]", "[g=0.75]", "[g=0.8]", "[g=0.85]", "[g=0.9]", "[g=0.95]")

perf_min <- min(sapply(X = res_gm, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gm, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gm, perf_min = perf_min-0.02, perf_max = perf_max+0.02)


############################################### GMP #################################################################

res_gmp <- list()

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_1_GMP", prob = 0.001, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_1 <- fit.classification(y=y, samples = samples, id = "result18_1_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_1, file=file.path('data/model/res_pa_GMP_18_1.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_1))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_2_GMP", prob = 0.001, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_2 <- fit.classification(y=y, samples = samples, id = "result18_2_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_2, file=file.path('data/model/res_pa_GMP_18_2.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_2))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_3_GMP", prob = 0.001, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_3 <- fit.classification(y=y, samples = samples, id = "result18_3_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_3, file=file.path('data/model/res_pa_GMP_18_3.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_3))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_4_GMP", prob = 0.001, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_4 <- fit.classification(y=y, samples = samples, id = "result18_4_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_4, file=file.path('data/model/res_pa_GMP_18_4.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_4))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_5_GMP", prob = 0.001, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_5 <- fit.classification(y=y, samples = samples, id = "result18_5_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_5, file=file.path('data/model/res_pa_GMP_18_5.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_5))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_6_GMP", prob = 0.01, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_6 <- fit.classification(y=y, samples = samples, id = "result18_6_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_6, file=file.path('data/model/res_pa_GMP_18_6.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_6))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_7_GMP", prob = 0.01, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_7 <- fit.classification(y=y, samples = samples, id = "result18_7_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_7, file=file.path('data/model/res_pa_GMP_18_7.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_7))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_8_GMP", prob = 0.01, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_8 <- fit.classification(y=y, samples = samples, id = "result18_8_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_8, file=file.path('data/model/res_pa_GMP_18_8.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_8))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_9_GMP", prob = 0.01, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_9 <- fit.classification(y=y, samples = samples, id = "result18_9_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_9, file=file.path('data/model/res_pa_GMP_18_9.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_9))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_10_GMP", prob = 0.01, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_10 <- fit.classification(y=y, samples = samples, id = "result18_10_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_10, file=file.path('data/model/res_pa_GMP_18_10.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_10))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_11_GMP", prob = 0.2, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_11 <- fit.classification(y=y, samples = samples, id = "result18_11_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_11, file=file.path('data/model/res_pa_GMP_18_11.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_11))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_12_GMP", prob = 0.2, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_12 <- fit.classification(y=y, samples = samples, id = "result18_12_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_12, file=file.path('data/model/res_pa_GMP_18_12.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_12))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_13_GMP", prob = 0.2, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_13 <- fit.classification(y=y, samples = samples, id = "result18_13_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_13, file=file.path('data/model/res_pa_GMP_18_13.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_13))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_14_GMP", prob = 0.2, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_14 <- fit.classification(y=y, samples = samples, id = "result18_14_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_14, file=file.path('data/model/res_pa_GMP_18_14.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_14))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_15_GMP", prob = 0.2, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_15 <- fit.classification(y=y, samples = samples, id = "result18_15_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_15, file=file.path('data/model/res_pa_GMP_18_15.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_15))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_16_GMP", prob = 0.4, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_16 <- fit.classification(y=y, samples = samples, id = "result18_16_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_16, file=file.path('data/model/res_pa_GMP_18_16.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_16))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_17_GMP", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_17 <- fit.classification(y=y, samples = samples, id = "result18_17_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_17, file=file.path('data/model/res_pa_GMP_18_17.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_17))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_18_GMP", prob = 0.4, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_18 <- fit.classification(y=y, samples = samples, id = "result18_18_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_18, file=file.path('data/model/res_pa_GMP_18_18.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_18))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_19_GMP", prob = 0.4, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_19 <- fit.classification(y=y, samples = samples, id = "result18_19_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_19, file=file.path('data/model/res_pa_GMP_18_19.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_19))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_20_GMP", prob = 0.4, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_20 <- fit.classification(y=y, samples = samples, id = "result18_20_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_20, file=file.path('data/model/res_pa_GMP_18_20.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_20))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_21_GMP", prob = 0.6, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_21 <- fit.classification(y=y, samples = samples, id = "result18_21_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_21, file=file.path('data/model/res_pa_GMP_18_21.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_21))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_22_GMP", prob = 0.6, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_22 <- fit.classification(y=y, samples = samples, id = "result18_22_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_22, file=file.path('data/model/res_pa_GMP_18_22.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_22))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_23_GMP", prob = 0.6, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_23 <- fit.classification(y=y, samples = samples, id = "result18_23_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_23, file=file.path('data/model/res_pa_GMP_18_23.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_23))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_24_GMP", prob = 0.6, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_24 <- fit.classification(y=y, samples = samples, id = "result18_24_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_24, file=file.path('data/model/res_pa_GMP_18_24.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_24))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_25_GMP", prob = 0.6, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_25 <- fit.classification(y=y, samples = samples, id = "result18_25_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_25, file=file.path('data/model/res_pa_GMP_18_25.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_25))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_26_GMP", prob = 0.8, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_26 <- fit.classification(y=y, samples = samples, id = "result18_26_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_26, file=file.path('data/model/res_pa_GMP_18_26.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_26))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_27_GMP", prob = 0.8, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_27 <- fit.classification(y=y, samples = samples, id = "result18_27_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_27, file=file.path('data/model/res_pa_GMP_18_27.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_27))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_28_GMP", prob = 0.8, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_28 <- fit.classification(y=y, samples = samples, id = "result18_28_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_28, file=file.path('data/model/res_pa_GMP_18_28.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_28))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_29_GMP", prob = 0.8, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_29 <- fit.classification(y=y, samples = samples, id = "result18_29_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_29, file=file.path('data/model/res_pa_GMP_18_29.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_29))


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_30_GMP", prob = 0.8, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_30 <- fit.classification(y=y, samples = samples, id = "result18_30_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 10)


save(res_pa_GMP_18_30, file=file.path('data/model/res_pa_GMP_18_30.RData'))

res_gmp <- c(res_gmp, list(res_pa_GMP_18_30))


#########################################################################################################################################
# plot

xlabs <- c("[p=0.001,g=0]", "[p=0.001,g=0.2]", "[p=0.001,g=0.4]", "[p=0.001,g=0.6]", "[p=0.001,g=0.8]",
           "[p=0.01,g=0]", "[p=0.01,g=0.2]", "[p=0.01,g=0.4]", "[p=0.01,g=0.6]", "[p=0.01,g=0.8]",
           "[p=0.2,g=0]", "[p=0.2,g=0.2]", "[p=0.2,g=0.4]", "[p=0.2,g=0.6]", "[p=0.2,g=0.8]", 
           "[p=0.4,g=0]", "[p=0.4,g=0.2]", "[p=0.4,g=0.4]", "[p=0.4,g=0.6]", "[p=0.4,g=0.8]", 
           "[p=0.6,g=0]", "[p=0.6,g=0.2]", "[p=0.6,g=0.4]", "[p=0.6,g=0.6]", "[p=0.6,g=0.8]", 
           "[p=0.8,g=0]", "[p=0.8,g=0.2]", "[p=0.8,g=0.4]", "[p=0.8,g=0.6]", "[p=0.8,g=0.8]")

# Plot for GMP models
title <- c("Result 18_GMP")
perf_min <- min(sapply(X = res_gmp, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmp, FUN = function(x){max(x$results$Accuracy)}))
perf_facet_boxplot(title, xlabs, res_gmp, perf_min = perf_min-0.01, perf_max = perf_max+0.01, perf_max)

