
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
library(doParallel)

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

# Plot for Result 25

for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_', id_list[i], '_LOOCV.RData'), collapse = ''))
}

for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GM_', id_list[i], '_LOOCV.RData'), collapse = ''))
}

res_model_25_5 <- list(res_pa_GM_25_5_LOOCV, res_pa_GMR_25_5_LOOCV)
res_model_25_10 <- list(res_pa_GM_25_10_LOOCV, res_pa_GMR_25_10_LOOCV)
res_model_25_15 <- list(res_pa_GM_25_15_LOOCV, res_pa_GMR_25_15_LOOCV)

title <- c("Result 25_5")
xlabs <- c("GM", "GMR")
perf_min <- min(sapply(X = res_model_25_5, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_5, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_5, perf_min = 0.65, perf_max = 0.71)


title <- c("Result 25_10")
xlabs <- c("GM", "GMR")
perf_min <- min(sapply(X = res_model_25_10, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_10, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_10, perf_min = 0.65, perf_max = 0.71)


title <- c("Result 25_15")
xlabs <- c("GM", "GMR")
perf_min <- min(sapply(X = res_model_25_15, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_15, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_15, perf_min = 0.65, perf_max = 0.71)


