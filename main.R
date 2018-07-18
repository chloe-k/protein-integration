
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




################################################### Result25  #################################################
res_models <- list()
#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result25_GMR", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_25 <- fit.classification(y=y, samples = samples, id = "result25_GMR", datapath = datapath, respath = respath,
                                    profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                    numTops=50)

save(res_pa_GMR_25, file=file.path('data/model/res_pa_GMR_25.RData'))

res_models <- c(res_models, list(res_pa_GMR_25))


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

# fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
#                id = "result25_GMR_d", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_25 <- fit.classification(y=y, samples = samples, id = "result25_GMR_d", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      numTops=50)

save(res_pa_GMR_d_25, file=file.path('data/model/res_pa_GMR_d_25.RData'))

res_models <- c(res_models, list(res_pa_GMR_d_25))



#########################################################################################################################################
# plot

# GM
title <- c("Result 25(15 intersect genes)")
xlabs <- c("GMR", "GMR_d")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.02, perf_max = perf_max+0.02)
