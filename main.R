
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

################################## Result 18_GMR ############################################################

num_cores <- detectCores()/2
registerDoParallel(cores = num_cores)


id_list <- c("18_1", "18_2", "18_3", "18_4", "18_5")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")
res_gmr_LOOCV <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i])
}

res_gmr_LOOCV <- list(res_pa_GMR_18_1_LOOCV, res_pa_GMR_18_2_LOOCV, res_pa_GMR_18_3_LOOCV, res_pa_GMR_18_4_LOOCV, res_pa_GMR_18_5_LOOCV)

title <- c("Result 18_GMR_LOOCV")
xlabs <- c("[g=0]", "[g=0.2]", "[g=0.4]", "[g=0.6]", "[g=0.8]")
perf_min <- min(sapply(X = res_gmr_LOOCV, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_LOOCV, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr_LOOCV, perf_min = perf_min-0.15, perf_max = perf_max+0.15)

#########################################################################################################################################
# Plot for GMR models

res_gmr <- list(res_pa_GMR_18_1, res_pa_GMR_18_2, res_pa_GMR_18_3, res_pa_GMR_18_4, res_pa_GMR_18_5)

title <- c("Result 18_GMR")
xlabs <- c("[g=0]", "[g=0.2]", "[g=0.4]", "[g=0.6]", "[g=0.8]")
perf_min <- min(sapply(X = res_gmr, FUN = function(x){mean(x$resample$Accuracy)}))
perf_max <- max(sapply(X = res_gmr, FUN = function(x){mean(x$resample$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr, perf_min = perf_min-0.15, perf_max = perf_max+0.15)

# Accuracy((A+D)/(A+B+C+D))
i=0
for(model in res_gmr){
  print(i)
  print(confusionMatrix(model, "none"))
  i <- i+1
}
