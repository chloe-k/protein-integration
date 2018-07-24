
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
################################## Result 18 in GMR ############################################################
num_cores <- detectCores()/2
registerDoParallel(cores = 4)

id_list <- c("18_1", "18_2", "18_3", "18_4", "18_5",
             "18_6", "18_7", "18_8")

Gamma_list <- c(0, 0.2, 0.4, 0.6, 0.8, 0.85, 0.9, 0.95)

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_18 <- foreach(i=6:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = "gmp", prob = 0.001, Gamma = Gamma_list[i])
}

for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_18_', i, '_LOOCV.RData'), collapse = ''))
}

# Plot for GMR models

res_gmr <- list(res_pa_GMR_18_1_LOOCV, res_pa_GMR_18_2_LOOCV, res_pa_GMR_18_3_LOOCV, res_pa_GMR_18_4_LOOCV, res_pa_GMR_18_5_LOOCV,
                res_pa_GMR_18_6_LOOCV, res_pa_GMR_18_7_LOOCV, res_pa_GMR_18_8_LOOCV)

title <- c("Result 18_GMR")
xlabs <- c("[g=0]", "[g=0.2]", "[g=0.4]", "[g=0.6]", "[g=0.8]", "[g=0.85]", "[g=0.9]", "[g=0.95]")
perf_min <- min(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr, perf_min = perf_min-0.15, perf_max = perf_max+0.15)

# ################################## Result 25 ############################################################
# 
# num_cores <- detectCores()/2
# registerDoParallel(cores = 4)
# 
# 
# id_list <- c("25_5", "25_10", "25_15")
# lim_list <- c(5, 10, 15)
# 
# 
# pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
#           "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")
# res_gmr_25 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
#   make_GMR_model(id=id_list[i], lim=lim_list[i], type_used = "gmp", prob = )
# }
# 
# res_gmr_d_25 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
#   make_GMR_d_model(id=id_list[i], lim = lim_list[i])
# }
# 
# 
# 
# 
# #########################################################################################################################################
# # Plot for Result 25
# 
# # res_models <- list(res_pa_GMR_25_5_LOOCV, res_pa_GMR_d_25_5_LOOCV, 
# #                    res_pa_GMR_25_10_LOOCV, res_pa_GMR_d_25_10_LOOCV,
# #                    res_pa_GMR_25_15_LOOCV, res_pa_GMR_d_25_15_LOOCV)
# # 
# # title <- c("Result 25")
# # xlabs <- c("GMR_5", "GMR_d_5", "GMR_10", "GMR_d_10", "GMR_15", "GMR_d_15")
# # perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
# # perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
# # perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.15, perf_max = perf_max+0.15)
# 
# for(i in 1:length(id_list)){
#   load(paste(c('data/model/res_pa_GMR_', id_list[i], '_LOOCV.RData'), collapse = ''))
# }
# 
# for(i in 1:length(id_list)){
#   load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '_LOOCV.RData'), collapse = ''))
# }
# 
# res_model_25_5 <- list(res_pa_GMR_25_5_LOOCV, res_pa_GMR_d_25_5_LOOCV)
# res_model_25_10 <- list(res_pa_GMR_25_10_LOOCV, res_pa_GMR_d_25_10_LOOCV)
# res_model_25_15 <- list(res_pa_GMR_25_15_LOOCV, res_pa_GMR_d_25_15_LOOCV)
# 
# title <- c("Result 25_5")
# xlabs <- c("GMR", "GMR_d")
# perf_min <- min(sapply(X = res_model_25_5, FUN = function(x){max(x$results$Accuracy)}))
# perf_max <- max(sapply(X = res_model_25_5, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_model_25_5, perf_min = perf_min-0.15, perf_max = perf_max+0.15)
# 
# 
# title <- c("Result 25_10")
# xlabs <- c("GMR", "GMR_d")
# perf_min <- min(sapply(X = res_model_25_10, FUN = function(x){max(x$results$Accuracy)}))
# perf_max <- max(sapply(X = res_model_25_10, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_model_25_10, perf_min = perf_min-0.15, perf_max = perf_max+0.15)
# 
# 
# title <- c("Result 25_15")
# xlabs <- c("GMR", "GMR_d")
# perf_min <- min(sapply(X = res_model_25_15, FUN = function(x){max(x$results$Accuracy)}))
# perf_max <- max(sapply(X = res_model_25_15, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_model_25_15, perf_min = perf_min-0.15, perf_max = perf_max+0.15)
