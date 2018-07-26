
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
perf_lineplot(title, xlabs, res_gmr, perf_min, perf_max, Gamma_list)

################################## Result 27 in GMR ############################################################
registerDoParallel(cores = 4)

id_list <- c("27_0.9_G", "27_0.9_M", "27_0.9_R", "27_0.9_GM", "27_0.9_GP", "27_0.9_MP", "27_0.9_GMP")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_27_0.9 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], prob = 0.001, Gamma = 0.9, type_used = type_list[i])
}




# write sigFeature
for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_', id_list[i], '.RData'), collapse = ''))
}

res_gmr_27 <- list(res_pa_GMR_27_0.9_G, res_pa_GMR_27_0.9_M, res_pa_GMR_27_0.9_R,
                       res_pa_GMR_27_0.9_GM, res_pa_GMR_27_0.9_GP, res_pa_GMR_27_0.9_MP,
                       res_pa_GMR_27_0.9_GMP)

profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR'), collapse = '')
  write.SigFeatures(res_fit=res_gmr_27[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

# Plot GMR
title <- c("Result 27 GMR")
xlabs <- rep(c("G", "M", "R", "GM", "GR", "MR", "GMR"), 5)
Gamma_list <- rep(c("0.2", "0.4", "0.6", "0.8", "0.9"), each=7)

id_list <- c("27_0.2_G", "27_0.2_M", "27_0.2_R", "27_0.2_GM", "27_0.2_GP", "27_0.2_MP", "27_0.2_GMP")
id_list <- c("27_0.4_G", "27_0.4_M", "27_0.4_R", "27_0.4_GM", "27_0.4_GP", "27_0.4_MP", "27_0.4_GMP")
id_list <- c("27_0.6_G", "27_0.6_M", "27_0.6_R", "27_0.6_GM", "27_0.6_GP", "27_0.6_MP", "27_0.6_GMP")
id_list <- c("27_0.8_G", "27_0.8_M", "27_0.8_R", "27_0.8_GM", "27_0.8_GP", "27_0.8_MP", "27_0.8_GMP")
id_list <- c("27_0.9_G", "27_0.9_M", "27_0.9_R", "27_0.9_GM", "27_0.9_GP", "27_0.9_MP", "27_0.9_GMP")

res_gmr <- list()
for(i in 1:length(id_list)){
  model<- get(load(paste(c('data/model/res_pa_GMR_', id_list[i], '.RData'), collapse = '')))
  res_gmr <- c(res_gmr, list(model))
}

perf_min <- min(sapply(X = res_gmr, FUN = function(x){mean(x$resample$Accuracy)}))
perf_max <- max(sapply(X = res_gmr, FUN = function(x){mean(x$resample$Accuracy)}))
# perf_boxplot(title, xlabs, res_gmr_27_0.2, perf_min = perf_min-0.01, perf_max = perf_max+0.01)
perf_lineplot(title, xlabs, res_gmr, perf_min, perf_max, Gamma_list)
