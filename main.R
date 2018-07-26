
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

################################## Result 25 ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = num_cores)


id_list <- c("25_10", "25_15")
lim_list <- c(10, 15)


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")
# res_gmr_25 <- foreach(i=1:1, .packages = pack) %dopar%{
#   make_GMR_model(id="25_10", lim=10)
# }

res_gmr_d_25 <- foreach(i=1:1, .packages = pack) %dopar%{
  make_GMR_d_model(id="25_15", lim = 15)
}




#########################################################################################################################################
# Plot for Result 25

res_models <- list(res_pa_GMR_25_5_LOOCV, res_pa_GMR_d_25_5_LOOCV, 
                   res_pa_GMR_25_10_LOOCV, res_pa_GMR_d_25_10_LOOCV,
                   res_pa_GMR_25_15_LOOCV, res_pa_GMR_d_25_15_LOOCV)

title <- c("Result 25")
xlabs <- c("GMR_5", "GMR_d_5", "GMR_10", "GMR_d_10", "GMR_15", "GMR_d_15")
perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.15, perf_max = perf_max+0.15)


res_model_25_5 <- list(res_pa_GMR_25_5_LOOCV, res_pa_GMR_d_25_5_LOOCV)
res_model_25_10 <- list(res_pa_GMR_25_10_LOOCV, res_pa_GMR_d_25_10_LOOCV)
res_model_25_15 <- list(res_pa_GMR_25_15_LOOCV, res_pa_GMR_d_25_15_LOOCV)

title <- c("Result 25")
xlabs <- c("GMR_5", "GMR_d_5")
perf_min <- min(sapply(X = res_model_25_5, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_5, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_5, perf_min = perf_min-0.15, perf_max = perf_max+0.15)


title <- c("Result 25_10")
xlabs <- c("GMR_10", "GMR_d_10")
perf_min <- min(sapply(X = res_model_25_10, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_10, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_10, perf_min = perf_min-0.15, perf_max = perf_max+0.15)


title <- c("Result 25_15")
xlabs <- c("GMR_15", "GMR_d_15")
perf_min <- min(sapply(X = res_model_25_15, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_15, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_15, perf_min = perf_min-0.15, perf_max = perf_max+0.15)
