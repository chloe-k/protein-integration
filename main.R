
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


# dppipath <- file.path(datapath, 'DppiGraph.rda')  # directed edge PPI
dppipath <- file.path(datapath, 'DppiGraph(Entrez).rda')  # directed edge PPI
#dppipath <- file.path(datapath, 'DppiGraph_rdc.rda')  # directed edge PPI
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

################################## Result 18_2 ############################################################

#------------------------- RNAseq(pathway) -------------------------#
testStatistic <- c("DESeq2")
profile_name <- c("rna(Entrez)")
x=list(rnaseq)

fit.iDRWPClass(x=x, y=y, globalGraph=g, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result18_2", 
               pranking = "t-test", mode = "G", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_G_18_2 <- fit.classification(y=y, samples = samples, id="result18_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                    method = "DRW", pranking = "t-test", classifier = "rf",
                                    nFolds = 5, numTops=50, iter = 50)

save(res_pa_G_18_2, file=file.path('data/model/res_pa_G_18_2.RData'))

summary(res_pa_G_18_2)
print(res_pa_G_18_2$results)
print(res_pa_G_18_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_G_18_2, id="result18_2", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- Methyl(pathway) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("meth(Entrez)")
x=list(imputed_methyl)


fit.iDRWPClass(x=x, y=y, globalGraph=m, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result18_2", 
               pranking = "t-test", mode = "M", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_M_18_2 <- fit.classification(y=y, samples = samples, id="result18_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                    method = "DRW", pranking = "t-test", classifier = "rf",
                                    nFolds = 5, numTops=50, iter = 50)

save(res_pa_M_18_2, file=file.path('data/model/res_pa_M_18_2.RData'))

summary(res_pa_M_18_2)
print(res_pa_M_18_2$results)
print(res_pa_M_18_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_M_18_2, id="result18_2", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using Pathway Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa(Entrez)")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=r, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result18_2", 
               pranking = "t-test", mode = "R", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_R_18_2 <- fit.classification(y=y, samples = samples, id="result18_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                    method = "DRW", pranking = "t-test", classifier = "rf",
                                    nFolds = 5, numTops=50, iter = 50)

save(res_pa_R_18_2, file=file.path('data/model/res_pa_R_18_2.RData'))

summary(res_pa_R_18_2)
print(res_pa_R_18_2$results)
print(res_pa_R_18_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_R_18_2, id="result18_2", profile_name=profile_name, method="DRW", respath=respath)

#------------------------- RPPA(pathway using PPI Graph) -------------------------#
testStatistic <- c("t-test")
profile_name <- c("rppa(Entrez)")
x=list(rppa)


fit.iDRWPClass(x=x, y=y, globalGraph=p, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result18_2", 
               pranking = "t-test", mode = "P", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_P_18_2 <- fit.classification(y=y, samples = samples, id="result18_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                    method = "DRW", pranking = "t-test", classifier = "rf",
                                    nFolds = 5, numTops=50, iter = 50)

save(res_pa_P_18_2, file=file.path('data/model/res_pa_P_18_2.RData'))

summary(res_pa_P_18_2)
print(res_pa_P_18_2$results)
print(res_pa_P_18_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_P_18_2, id="result18_2", profile_name=profile_name, method="DRW", respath=respath)

# plot
title <- c("Result 18_2")
xlabs <- c("G", "M", "R", "P")
res_models <- list(res_pa_G_18_2, res_pa_M_18_2, res_pa_R_18_2, res_pa_P_18_2)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_G_18_2$results$Accuracy[1])

