
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

res_gmr <- list()

#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_1_GMR", prob = 0.001, Gamma = 0, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_1 <- fit.classification(y=y, samples = samples, id = "result18_1_GMR", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_18_1, file=file.path('data/model/res_pa_GMR_18_1.RData'))

res_gmr <- c(res_gmr, list(res_pa_GMR_18_1))


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_2_GMR", prob = 0.001, Gamma = 0.2, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_2 <- fit.classification(y=y, samples = samples, id = "result18_2_GMR", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_18_2, file=file.path('data/model/res_pa_GMR_18_2.RData'))

res_gmr <- c(res_gmr, list(res_pa_GMR_18_2))


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_3_GMR", prob = 0.001, Gamma = 0.4, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_3 <- fit.classification(y=y, samples = samples, id = "result18_3_GMR", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_18_3, file=file.path('data/model/res_pa_GMR_18_3.RData'))

res_gmr <- c(res_gmr, list(res_pa_GMR_18_3))


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_4_GMR", prob = 0.001, Gamma = 0.6, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_4 <- fit.classification(y=y, samples = samples, id = "result18_4_GMR", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_18_4, file=file.path('data/model/res_pa_GMR_18_4.RData'))

res_gmr <- c(res_gmr, list(res_pa_GMR_18_4))


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_5_GMR", prob = 0.001, Gamma = 0.8, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_5 <- fit.classification(y=y, samples = samples, id = "result18_5_GMR", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_18_5, file=file.path('data/model/res_pa_GMR_18_5.RData'))

res_gmr <- c(res_gmr, list(res_pa_GMR_18_5))



############################################################################################
# plot

xlabs <- c("[g=0]", "[g=0.2]", "[g=0.4]", "[g=0.6]", "[g=0.8]")

# Plot for GMR models
title <- c("Result 18_GMR")
perf_min <- min(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr, perf_min = perf_min-0.02, perf_max = perf_max+0.02)


############################################### GMR_d ######################################################################

res_gmr_d <- list()

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_1_GMR_d", prob = 0.001, Gamma = 0, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_1 <- fit.classification(y=y, samples = samples, id = "result18_1_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_1, file=file.path('data/model/res_pa_GMR_d_18_1.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_1))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_2_GMR_d", prob = 0.001, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_2 <- fit.classification(y=y, samples = samples, id = "result18_2_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_2, file=file.path('data/model/res_pa_GMR_d_18_2.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_2))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_3_GMR_d", prob = 0.001, Gamma = 0.4, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_3 <- fit.classification(y=y, samples = samples, id = "result18_3_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_3, file=file.path('data/model/res_pa_GMR_d_18_3.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_3))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_4_GMR_d", prob = 0.001, Gamma = 0.6, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_4 <- fit.classification(y=y, samples = samples, id = "result18_4_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_4, file=file.path('data/model/res_pa_GMR_d_18_4.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_4))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_5_GMR_d", prob = 0.001, Gamma = 0.8, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_5 <- fit.classification(y=y, samples = samples, id = "result18_5_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_5, file=file.path('data/model/res_pa_GMR_d_18_5.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_5))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_6_GMR_d", prob = 0.01, Gamma = 0, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_6 <- fit.classification(y=y, samples = samples, id = "result18_6_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_6, file=file.path('data/model/res_pa_GMR_d_18_6.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_6))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_7_GMR_d", prob = 0.01, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_7 <- fit.classification(y=y, samples = samples, id = "result18_7_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_7, file=file.path('data/model/res_pa_GMR_d_18_7.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_7))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_8_GMR_d", prob = 0.01, Gamma = 0.4, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_8 <- fit.classification(y=y, samples = samples, id = "result18_8_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_8, file=file.path('data/model/res_pa_GMR_d_18_8.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_8))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_9_GMR_d", prob = 0.01, Gamma = 0.6, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_9 <- fit.classification(y=y, samples = samples, id = "result18_9_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_9, file=file.path('data/model/res_pa_GMR_d_18_9.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_9))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_10_GMR_d", prob = 0.01, Gamma = 0.8, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_10 <- fit.classification(y=y, samples = samples, id = "result18_10_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_10, file=file.path('data/model/res_pa_GMR_d_18_10.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_10))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_11_GMR_d", prob = 0.2, Gamma = 0, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_11 <- fit.classification(y=y, samples = samples, id = "result18_11_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_11, file=file.path('data/model/res_pa_GMR_d_18_11.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_11))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_12_GMR_d", prob = 0.2, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_12 <- fit.classification(y=y, samples = samples, id = "result18_12_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_12, file=file.path('data/model/res_pa_GMR_d_18_12.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_12))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_13_GMR_d", prob = 0.2, Gamma = 0.4, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_13 <- fit.classification(y=y, samples = samples, id = "result18_13_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_13, file=file.path('data/model/res_pa_GMR_d_18_13.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_13))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_14_GMR_d", prob = 0.2, Gamma = 0.6, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_14 <- fit.classification(y=y, samples = samples, id = "result18_14_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_14, file=file.path('data/model/res_pa_GMR_d_18_14.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_14))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_15_GMR_d", prob = 0.2, Gamma = 0.8, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_15 <- fit.classification(y=y, samples = samples, id = "result18_15_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_15, file=file.path('data/model/res_pa_GMR_d_18_15.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_15))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_16_GMR_d", prob = 0.4, Gamma = 0, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_16 <- fit.classification(y=y, samples = samples, id = "result18_16_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_16, file=file.path('data/model/res_pa_GMR_d_18_16.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_16))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_17_GMR_d", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_17 <- fit.classification(y=y, samples = samples, id = "result18_17_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_17, file=file.path('data/model/res_pa_GMR_d_18_17.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_17))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_18_GMR_d", prob = 0.4, Gamma = 0.4, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_18 <- fit.classification(y=y, samples = samples, id = "result18_18_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_18, file=file.path('data/model/res_pa_GMR_d_18_18.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_18))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_19_GMR_d", prob = 0.4, Gamma = 0.6, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_19 <- fit.classification(y=y, samples = samples, id = "result18_19_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_19, file=file.path('data/model/res_pa_GMR_d_18_19.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_19))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_20_GMR_d", prob = 0.4, Gamma = 0.8, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_20 <- fit.classification(y=y, samples = samples, id = "result18_20_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_20, file=file.path('data/model/res_pa_GMR_d_18_20.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_20))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_21_GMR_d", prob = 0.6, Gamma = 0, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_21 <- fit.classification(y=y, samples = samples, id = "result18_21_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_21, file=file.path('data/model/res_pa_GMR_d_18_21.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_21))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_22_GMR_d", prob = 0.6, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_22 <- fit.classification(y=y, samples = samples, id = "result18_22_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_22, file=file.path('data/model/res_pa_GMR_d_18_22.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_22))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_23_GMR_d", prob = 0.6, Gamma = 0.4, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_23 <- fit.classification(y=y, samples = samples, id = "result18_23_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_23, file=file.path('data/model/res_pa_GMR_d_18_23.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_23))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_24_GMR_d", prob = 0.6, Gamma = 0.6, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_24 <- fit.classification(y=y, samples = samples, id = "result18_24_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_24, file=file.path('data/model/res_pa_GMR_d_18_24.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_24))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_25_GMR_d", prob = 0.6, Gamma = 0.8, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_25 <- fit.classification(y=y, samples = samples, id = "result18_25_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_25, file=file.path('data/model/res_pa_GMR_d_18_25.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_25))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_26_GMR_d", prob = 0.8, Gamma = 0, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_26 <- fit.classification(y=y, samples = samples, id = "result18_26_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_26, file=file.path('data/model/res_pa_GMR_d_18_26.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_26))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_27_GMR_d", prob = 0.8, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_27 <- fit.classification(y=y, samples = samples, id = "result18_27_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_27, file=file.path('data/model/res_pa_GMR_d_18_27.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_27))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_28_GMR_d", prob = 0.8, Gamma = 0.4, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_28 <- fit.classification(y=y, samples = samples, id = "result18_28_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_28, file=file.path('data/model/res_pa_GMR_d_18_28.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_28))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_29_GMR_d", prob = 0.8, Gamma = 0.6, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_29 <- fit.classification(y=y, samples = samples, id = "result18_29_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_29, file=file.path('data/model/res_pa_GMR_d_18_29.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_29))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_30_GMR_d", prob = 0.8, Gamma = 0.8, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_30 <- fit.classification(y=y, samples = samples, id = "result18_30_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_d_18_30, file=file.path('data/model/res_pa_GMR_d_18_30.RData'))

res_gmr_d <- c(res_gmr_d, list(res_pa_GMR_d_18_30))

#########################################################################################################################################
# plot

xlabs <- c("[p=0.001,g=0]", "[p=0.001,g=0.2]", "[p=0.001,g=0.4]", "[p=0.001,g=0.6]", "[p=0.001,g=0.8]",
           "[p=0.01,g=0]", "[p=0.01,g=0.2]", "[p=0.01,g=0.4]", "[p=0.01,g=0.6]", "[p=0.01,g=0.8]",
           "[p=0.2,g=0]", "[p=0.2,g=0.2]", "[p=0.2,g=0.4]", "[p=0.2,g=0.6]", "[p=0.2,g=0.8]", 
           "[p=0.4,g=0]", "[p=0.4,g=0.2]", "[p=0.4,g=0.4]", "[p=0.4,g=0.6]", "[p=0.4,g=0.8]", 
           "[p=0.6,g=0]", "[p=0.6,g=0.2]", "[p=0.6,g=0.4]", "[p=0.6,g=0.6]", "[p=0.6,g=0.8]", 
           "[p=0.8,g=0]", "[p=0.8,g=0.2]", "[p=0.8,g=0.4]", "[p=0.8,g=0.6]", "[p=0.8,g=0.8]")

# Plot for GMR_d models
title <- c("Result 18_GMR_d")
perf_min <- min(sapply(X = res_gmr_d, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_d, FUN = function(x){max(x$results$Accuracy)}))
perf_facet_boxplot(title, xlabs, res_gmr_d, perf_min = perf_min-0.01, perf_max = perf_max+0.01, perf_max)
