
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
# ################################## Result 23 in GMR_d ############################################################
# num_cores <- detectCores()/2
# registerDoParallel(cores = num_cores)
# 
# id_list <- c("23_G", "23_M", "23_P", "23_GM", "23_GP", "23_MP", "23_GMP")
# type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")
# 
# 
# pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
#           "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")
# 
# res_gmr_d_23 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
#   make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.4, Gamma = 0.2)
# }


# ################################## Result 18 in GMR ############################################################
# num_cores <- detectCores()/2
# registerDoParallel(cores = 4)
# 
# id_list <- c("18_1", "18_2", "18_3", "18_4", "18_5")
# 
# Gamma_list <- c(0, 0.2, 0.4, 0.6, 0.8)
# # make_GMR_model(id=id_list[1], type_used = "gmp", prob = 0.001, Gamma = Gamma_list[1])
# 
# 
# # make_GMR_model(id="18_1", type_used = "gmp", prob = 0.001, Gamma = 0)
# 
# pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
#           "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")
# 
# res_gmr_18 <- foreach(i=2:length(id_list), .packages = pack) %dopar%{
#   make_GMR_model(id=id_list[i], type_used = "gmp", prob = 0.001, Gamma = Gamma_list[i])
# }
# 
# for(i in 1:length(id_list)){
#   load(paste(c('data/model/res_pa_GMR_18_', i, '_LOOCV.RData'), collapse = ''))
# }
# 
# 
# # Plot for GMR models
# 
# res_gmr <- list(res_pa_GMR_18_1_LOOCV, res_pa_GMR_18_2_LOOCV, res_pa_GMR_18_3_LOOCV, 
#                 res_pa_GMR_18_4_LOOCV, res_pa_GMR_18_5_LOOCV)
# 
# title <- c("Result 18_GMR")
# xlabs <- c("[g=0]", "[g=0.2]", "[g=0.4]", "[g=0.6]", "[g=0.8]")
# perf_min <- min(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
# perf_max <- max(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_gmr, perf_min = perf_min-0.15, perf_max = perf_max+0.15)
# 
# 
# #########################################################################################################################################
# # plot
# 
# # GM
# res_gm_23 <- list(res_pa_GM_23_G, res_pa_GM_23_M)
# 
# title <- c("Result 23_GM")
# xlabs <- c("G", "M")
# 
# perf_min <- min(sapply(X = res_gm_23, FUN = function(x){max(x$results$Accuracy)}))
# perf_max <- max(sapply(X = res_gm_23, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_gm_23, perf_min = perf_min-0.02, perf_max = perf_max+0.02)
# 
# #GMR_d
# res_gmr_d_23 <- list(res_pa_GMR_d_23_G_LOOCV, res_pa_GMR_d_23_M_LOOCV, res_pa_GMR_d_23_P_LOOCV, 
#                      res_pa_GMR_d_23_GM_LOOCV, res_pa_GMR_d_23_GP_LOOCV, res_pa_GMR_d_23_MP_LOOCV, 
#                      res_pa_GMR_d_23_GMP_LOOCV)
# 
# title <- c("Result 23_GMR_d")
# xlabs <- c("G", "M", "P", "GM", "GP", "MP", "GMP")
# 
# perf_min <- min(sapply(X = res_gmr_d_23, FUN = function(x){max(x$results$Accuracy)}))
# perf_max <- max(sapply(X = res_gmr_d_23, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_gmr_d_23, perf_min = perf_min-0.02, perf_max = perf_max+0.02)
# 

################################## Result 18_GMR_d ############################################################
num_cores <- detectCores()/2
registerDoParallel(cores = 4)


id_list <- c("18_1", "18_2", "18_3", "18_4", "18_5",
             "18_6", "18_7", "18_8", "18_9", "18_10",
             "18_11", "18_12", "18_13", "18_14", "18_15",
             "18_16", "18_17", "18_18", "18_19", "18_20",
             "18_21", "18_22", "18_23", "18_24", "18_25",
             "18_26", "18_27", "18_28", "18_29", "18_30")


prob_list <- c(0.001, 0.001, 0.001, 0.001, 0.001,
               0.01, 0.01, 0.01, 0.01, 0.01,
               0.2, 0.2, 0.2, 0.2, 0.2,
               0.4, 0.4, 0.4, 0.4, 0.4,
               0.6, 0.6, 0.6, 0.6, 0.6,
               0.8, 0.8, 0.8, 0.8, 0.8)

Gamma_list <- c(0, 0.2, 0.4, 0.6, 0.8,
                0, 0.2, 0.4, 0.6, 0.8,
                0, 0.2, 0.4, 0.6, 0.8,
                0, 0.2, 0.4, 0.6, 0.8,
                0, 0.2, 0.4, 0.6, 0.8,
                0, 0.2, 0.4, 0.6, 0.8)

make_GMR_d_model(id=id_list[1], type_used = "gmp", prob = prob_list[1], Gamma = Gamma_list[1])

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_LOOCV <- foreach(i=2:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = "gmp", prob = prob_list[i], Gamma = Gamma_list[i])
}

for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_d_18_', i, '_LOOCV.RData'), collapse = ''))
}

##########################################################################################
# Plot for GMR_d models



res_gmr_d_LOOCV <- list(res_pa_GMR_d_18_1_LOOCV, res_pa_GMR_d_18_2_LOOCV, res_pa_GMR_d_18_3_LOOCV, res_pa_GMR_d_18_4_LOOCV, res_pa_GMR_d_18_5_LOOCV,
                  res_pa_GMR_d_18_6_LOOCV, res_pa_GMR_d_18_7_LOOCV, res_pa_GMR_d_18_8_LOOCV, res_pa_GMR_d_18_9_LOOCV, res_pa_GMR_d_18_10_LOOCV,
                  res_pa_GMR_d_18_11_LOOCV, res_pa_GMR_d_18_12_LOOCV, res_pa_GMR_d_18_13_LOOCV, res_pa_GMR_d_18_14_LOOCV, res_pa_GMR_d_18_15_LOOCV,
                  res_pa_GMR_d_18_16_LOOCV, res_pa_GMR_d_18_17_LOOCV, res_pa_GMR_d_18_18_LOOCV, res_pa_GMR_d_18_19_LOOCV, res_pa_GMR_d_18_20_LOOCV,
                  res_pa_GMR_d_18_21_LOOCV, res_pa_GMR_d_18_22_LOOCV, res_pa_GMR_d_18_23_LOOCV, res_pa_GMR_d_18_24_LOOCV, res_pa_GMR_d_18_25_LOOCV,
                  res_pa_GMR_d_18_26_LOOCV, res_pa_GMR_d_18_27_LOOCV, res_pa_GMR_d_18_28_LOOCV, res_pa_GMR_d_18_29_LOOCV, res_pa_GMR_d_18_30_LOOCV)

title <- c("Result 18_GMR_d_LOOCV")
xlabs <- c("[p=0.001,g=0]", "[p=0.001,g=0.2]", "[p=0.001,g=0.4]", "[p=0.001,g=0.6]", "[p=0.001,g=0.8]",
           "[p=0.01,g=0]", "[p=0.01,g=0.2]", "[p=0.01,g=0.4]", "[p=0.01,g=0.6]", "[p=0.01,g=0.8]",
           "[p=0.2,g=0]", "[p=0.2,g=0.2]", "[p=0.2,g=0.4]", "[p=0.2,g=0.6]", "[p=0.2,g=0.8]",
           "[p=0.4,g=0]", "[p=0.4,g=0.2]", "[p=0.4,g=0.4]", "[p=0.4,g=0.6]", "[p=0.4,g=0.8]",
           "[p=0.6,g=0]", "[p=0.6,g=0.2]", "[p=0.6,g=0.4]", "[p=0.6,g=0.6]", "[p=0.6,g=0.8]",
           "[p=0.8,g=0]", "[p=0.8,g=0.2]", "[p=0.8,g=0.4]", "[p=0.8,g=0.6]", "[p=0.8,g=0.8]")

perf_min <- min(sapply(X = res_gmr_d_LOOCV, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_d_LOOCV, FUN = function(x){max(x$results$Accuracy)}))
perf_facet_boxplot(title, xlabs, res_gmr_d_LOOCV, perf_min = perf_min-0.05, perf_max = perf_max+0.05, perf_max)

