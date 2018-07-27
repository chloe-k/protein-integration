
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

# p <- DppiGraph
# V(p)$name <-paste("p",V(p)$name,sep="")

y=list(good_samples, poor_samples)
#---------------DRW---------------#
# RNAseq pathway profile


# Methylation pathway profile


# RPPA pathway profile


#----------------------------------------iDRW-----------------------------------------------------------#

################################## Result 29 in GR_d ############################################################
registerDoParallel(cores = 8)
id_list <- c("29_1", "29_2", "29_3", "29_4",
             "29_5", "29_6", "29_7", "29_8",
             "29_9", "29_10", "29_11", "29_12",
             "29_13", "29_14", "29_15", "29_16")


id_list <- c("28_0.2", "28_0.4", "28_0.6", "28_0.8", "28_0.9")
Gamma_list <- c(0.2, 0.4, 0.6, 0.8, 0.9)
make_GR_model(id=id_list[1], prob = 0.001, Gamma = Gamma_list[1])
pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gr_28 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GR_model(id=id_list[i], prob = 0.001, Gamma = Gamma_list[i])
}


res_models <- list()
for(i in 1:length(id_list)){

  load(paste(c('data/model/res_pa_GR_', id_list[i], '_LOOCV.RData'), collapse = ''))
}

res_gr_28 <- list(res_pa_GR_28_0.2_LOOCV, res_pa_GR_28_0.4_LOOCV, res_pa_GR_28_0.6_LOOCV,
                  res_pa_GR_28_0.8_LOOCV, res_pa_GR_28_0.9_LOOCV)

# profile_name <- c("rna(Entrez)", "rppa(Entrez)")
# for(i in 1:length(id_list)){
#   result_name <- paste(c('result',id_list[i],'_GR'), collapse = '')
#   write.SigFeatures(res_fit=res_gr_28[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
# }

#########################################################################
# Plot GR_d
title <- c("Result 28 GR_d")
xlabs <- c("[p=0.2,g=0.2]", "[p=0.2,g=0.4]", "[p=0.2,g=0.6]", "[p=0.2,g=0.8]",
           "[p=0.4,g=0.2]", "[p=0.4,g=0.4]", "[p=0.4,g=0.6]", "[p=0.4,g=0.8]",
           "[p=0.6,g=0.2]", "[p=0.6,g=0.4]", "[p=0.6,g=0.6]", "[p=0.6,g=0.8]",
           "[p=0.8,g=0.2]", "[p=0.8,g=0.4]", "[p=0.8,g=0.6]", "[p=0.8,g=0.8]")
prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), each=4)


perf_min <- min(sapply(X = res_gr_28, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gr_28, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gr_28, perf_min = perf_min-0.05, perf_max = perf_max+0.05)
# perf_lineplot(title, xlabs, res_gr_28, perf_min=55, perf_max=95)
perf_facet_boxplot(title, xlabs, res_models, perf_min = perf_min-0.15, perf_max = perf_max+0.15, perf_max, prob_list = prob_list, Gamma_list = Gamma_list)
perf_heatmap(title, xlabs, res_models)
