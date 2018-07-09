
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


<<<<<<< HEAD
# dppipath <- file.path(datapath, 'DppiGraph(Entrez).rda')  # directed edge PPI
# dppipath <- file.path(datapath, 'DppiGraph_rdc.rda')  # directed edge PPI
dppipath <- file.path(datapath, 'DppiGraph_W_str.rda')

=======
# dppipath <- file.path(datapath, 'DppiGraph.rda')  # directed edge PPI
dppipath <- file.path(datapath, 'DppiGraph(Entrez).rda')  # directed edge PPI
#dppipath <- file.path(datapath, 'DppiGraph_rdc.rda')  # directed edge PPI
>>>>>>> 25cb2006a4f48f4aeae0648fe7469fe5df5600b3
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


################################## Result 22 ############################################################

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
dppipath <- file.path(datapath, 'DppiGraph_W_str.rda')  # directed edge PPI
load(file.path(dppipath))

p <- DppiGraph
V(p)$name <-paste("p",V(p)$name,sep="")

testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result22",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_22 <- fit.classification(y=y, samples = samples, id = "result22", datapath = datapath, respath = respath, profile_name = profile_name,
                                    method = "DRW", pranking = "t-test", classifier = "rf",
                                    nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_22, file=file.path('data/model/res_pa_GMP_22.RData'))

summary(res_pa_GMP_22)
print(res_pa_GMP_22$results)

write.SigFeatures(res_fit=res_pa_GMP_22, id = "result22", profile_name=profile_name, method="DRW", respath=respath)


# plot
title <- c("Result 22")
xlabs <- c("GM", "GMR", "GMR_d", "GMP(unweighted STRING PPI)", "GMP(weighted STRING PPI)")
res_models <- list(res_pa_GM_RF_13, res_pa_GMR_RF_13, res_pa_GMR_d_RF_13, res_pa_GMP_21, res_pa_GMP_22)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_RF_13$results$Accuracy[1])
