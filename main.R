
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

################################## Result 18_GMR_d ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = num_cores)


id_list <- c("18_1", "18_2", "18_3", "18_4", "18_5",
             "18_6", "18_7", "18_8", "18_9", "18_10",
             "18_11", "18_12", "18_13", "18_14", "18_15",
             "18_16", "18_17", "18_18", "18_19", "18_20",
             "18_21", "18_22", "18_23", "18_24", "18_25",
             "18_26", "18_27", "18_28", "18_29", "18_30")

# make_GMP_model('18_1')
pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")
res_gmr_d <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i])
}


res_models <- list(res_pa_GMR_d_18_1, res_pa_GMR_d_18_2, res_pa_GMR_d_18_3, res_pa_GMR_d_18_4, res_pa_GMR_d_18_5,
                  res_pa_GMR_d_18_6, res_pa_GMR_d_18_7, res_pa_GMR_d_18_8, res_pa_GMR_d_18_9, res_pa_GMR_d_18_10,
                  res_pa_GMR_d_18_11, res_pa_GMR_d_18_12, res_pa_GMR_d_18_13, res_pa_GMR_d_18_14, res_pa_GMR_d_18_15,
                  res_pa_GMR_d_18_16, res_pa_GMR_d_18_17, res_pa_GMR_d_18_18, res_pa_GMR_d_18_19, res_pa_GMR_d_18_20,
                  res_pa_GMR_d_18_21, res_pa_GMR_d_18_22, res_pa_GMR_d_18_23, res_pa_GMR_d_18_24, res_pa_GMR_d_18_25,
                  res_pa_GMR_d_18_26, res_pa_GMR_d_18_27, res_pa_GMR_d_18_28, res_pa_GMR_d_18_29, res_pa_GMR_d_18_30)
#########################################################################################################################################
# plot
# 
# xlabs <- c("[p=0.001,g=0]", "[p=0.001,g=0.2]", "[p=0.001,g=0.4]", "[p=0.001,g=0.6]", "[p=0.001,g=0.8]",
#            "[p=0.01,g=0]", "[p=0.01,g=0.2]", "[p=0.01,g=0.4]", "[p=0.01,g=0.6]", "[p=0.01,g=0.8]",
#            "[p=0.2,g=0]", "[p=0.2,g=0.2]", "[p=0.2,g=0.4]", "[p=0.2,g=0.6]", "[p=0.2,g=0.8]",
#            "[p=0.4,g=0]", "[p=0.4,g=0.2]", "[p=0.4,g=0.4]", "[p=0.4,g=0.6]", "[p=0.4,g=0.8]",
#            "[p=0.6,g=0]", "[p=0.6,g=0.2]", "[p=0.6,g=0.4]", "[p=0.6,g=0.6]", "[p=0.6,g=0.8]",
#            "[p=0.8,g=0]", "[p=0.8,g=0.2]", "[p=0.8,g=0.4]", "[p=0.8,g=0.6]", "[p=0.8,g=0.8]")
# 
# # Plot for GMP models
# title <- c("Result 18_GMR_d")
# perf_min <- min(sapply(X = res_models, FUN = function(x){mean(x$resample$Accuracy)}))
# perf_max <- max(sapply(X = res_models, FUN = function(x){mean(x$resample$Accuracy)}))
# perf_facet_boxplot(title, xlabs, res_models, perf_min = perf_min-0.15, perf_max = perf_max+0.15, perf_max)
