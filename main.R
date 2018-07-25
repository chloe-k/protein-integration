
########################### package load ############################
# source("https://bioconductor.org/biocLite.R")

library(KEGGgraph)
library(igraph)
library(ggplot2)
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
title <- c("Result 18 GMR_d")
perf_heatmap(title, xlabs, res_gmr_d)

################################## Result 23 in GMR_d ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = 4)

id_list <- c("23_1", "23_2", "23_3", "23_4", "23_5", "23_6", "23_7")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


make_GMR_d_model(id=id_list[1], type_used = type_list[1], prob = 0.4, Gamma = 0.2)

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_23 <- foreach(i=2:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.4, Gamma = 0.2)
}


for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_d_23_', i, '_LOOCV.RData'), collapse = ''))
}

res_gmr_d <- list(res_pa_GMR_d_23_1_LOOCV, res_pa_GMR_d_23_2_LOOCV, res_pa_GMR_d_23_3_LOOCV, res_pa_GMR_d_23_4_LOOCV, 
                  res_pa_GMR_d_23_5_LOOCV, res_pa_GMR_d_23_6_LOOCV, res_pa_GMR_d_23_7_LOOCV)

# Plot for GMR_d model
title <- c("Result 23_GMR_d")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr_d, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_d, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr_d, perf_min = perf_min-0.005, perf_max = perf_max+0.005)
 
