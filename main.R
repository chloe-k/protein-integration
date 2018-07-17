
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




################################################### Result18_all  #################################################
# #------------------------- RNAseq + Methyl -------------------------#
# gm <- g %du% m
# testStatistic <- c("DESeq2", "t-test")
# profile_name <- c("rna(Entrez)", "meth(Entrez)")
# x=list(rnaseq, imputed_methyl)
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_0.75_GM", prob = 0.001, Gamma = 0.75, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)
# 
# res_pa_GM_18_0.75 <- fit.classification(y=y, samples = samples, id = "result18_0.75_GM", datapath = datapath, respath = respath,
#                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                                      nFolds = 5, numTops=50, iter = 20)
# 
# save(res_pa_GM_18_0.75, file=file.path('data/model/res_pa_GM_18_0.75.RData'))
# 
# #------------------------- RNAseq + Methyl -------------------------#
# gm <- g %du% m
# testStatistic <- c("DESeq2", "t-test")
# profile_name <- c("rna(Entrez)", "meth(Entrez)")
# x=list(rnaseq, imputed_methyl)
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_0.9_GM", prob = 0.001, Gamma = 0.9, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)
# 
# res_pa_GM_18_0.9 <- fit.classification(y=y, samples = samples, id = "result18_0.9_GM", datapath = datapath, respath = respath,
#                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                                         nFolds = 5, numTops=50, iter = 20)
# 
# save(res_pa_GM_18_0.9, file=file.path('data/model/res_pa_GM_18_0.9.RData'))
# 
# 
# #------------------------- RNAseq + Methyl -------------------------#
# gm <- g %du% m
# testStatistic <- c("DESeq2", "t-test")
# profile_name <- c("rna(Entrez)", "meth(Entrez)")
# x=list(rnaseq, imputed_methyl)
# 
# fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result18_0.95_GM", prob = 0.001, Gamma = 0.95, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)
# 
# res_pa_GM_18_0.95 <- fit.classification(y=y, samples = samples, id = "result18_0.95_GM", datapath = datapath, respath = respath,
#                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                                         nFolds = 5, numTops=50, iter = 20)
# 
# save(res_pa_GM_18_0.95, file=file.path('data/model/res_pa_GM_18_0.95.RData'))
# 
# res_pa_GMP_18_1 <- fit.classification(y=y, samples = samples, id = "result18_1_GMP", datapath = datapath, respath = respath,
#                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                                       nFolds = 5, numTops=50, iter = 20)


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
perf_min <- min(sapply(X = res_gmr_d, FUN = function(x){mean(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_d, FUN = function(x){mean(x$results$Accuracy)}))
perf_facet_boxplot(title, xlabs, res_gmr_d, perf_min = perf_min-0.01, perf_max = perf_max+0.01, perf_max)

# Plot for GMP models
title <- c("Result 18_GMP")
perf_min <- min(sapply(X = res_gmp, FUN = function(x){mean(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmp, FUN = function(x){mean(x$results$Accuracy)}))
perf_facet_boxplot(title, xlabs, res_gmp, perf_min = perf_min-0.02, perf_max = perf_max+0.02)