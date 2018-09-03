
########################### package load ############################
# source("https://bioconductor.org/biocLite.R")

library(KEGGgraph)
library(igraph)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(viridis)
library(reshape2)
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
library(doParallel)
library(ROCR)

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
# data_all_path <- file.path(datapath, "Dup_rppa_data.RData")
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



# dppipath <- file.path(datapath, 'DppiGraph(Entrez).rda')  # directed edge PPI
# # dppipath <- file.path(datapath, 'DppiGraph_rdc.rda')  # directed edge PPI
# # dppipath <- file.path(datapath, 'DppiGraph_W_str.rda')
# 
# if(!file.exists(dppipath)) {
#   print('ppiGraph does not exist, now creating ppiGraph start')
#   cons_ppi(datapath, gdacpath, rppa)
# }
# load(file.path(dppipath))


# directed pathway graph provided in DRWPClass
g <- directGraph 
V(g)$name <- paste("g",V(g)$name,sep="")

m <- directGraph
V(m)$name <-paste("m",V(m)$name,sep="")

r <- directGraph
V(r)$name <-paste("p",V(r)$name,sep="")

# p <- DppiGraph
# V(p)$name <-paste("p",V(p)$name,sep="")

y=list(good_samples, poor_samples)
#---------------DRW---------------#
# RNAseq pathway profile


# Methylation pathway profile


# RPPA pathway profile



#----------------------------------------iDRW-----------------------------------------------------------#
# GMP
# result_name <- paste(c('result','18_3','_GMR'), collapse = '')
# profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
# fit.classification(y=y, samples = samples, id = result_name, datapath = datapath, respath = respath,
#                    profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                    nFolds = 5, numTops=50, iter = 10)


# GP
result_name <- paste(c('result','28_0.6','_GR'), collapse = '')
profile_name <- c("rna(Entrez)", "rppa(Entrez)")
fit.classification(y=y, samples = samples, id = result_name, datapath = datapath, respath = respath,
                            profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                            nFolds = 5, numTops=50, iter = 10)