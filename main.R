
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

######################################################################

sapply(file.path("utils",list.files("utils", pattern="*.R")),source)

# make directory
datapath <- file.path('data')
if(!dir.exists(datapath)) dir.create(datapath)

gdacpath <- file.path(datapath, 'BRCA_GDAC')

respath <- file.path('result')
if(!dir.exists(respath)) dir.create(respath)


# read RNAseq, Methylation data, RPPA data
data_all_path <- file.path(datapath, "data.RData")
if(!file.exists(data_all_path)) {
  year <- 3
  read_data(year, datapath)
}
load(data_all_path)


# read rda data
graphpath <- file.path(datapath,'directGraph.rda')
pathSetpath <- file.path(datapath,'pathSet.rda')
if(!(file.exists(graphpath) && file.exists(pathSetpath))) {
  print('directGraph and pathSet do not exist, now creating directGraph and pathSet is start')
  cons_KEGGgraph(datapath)
}
load(file.path(graphpath))
load(file.path(pathSetpath))


ppipath <- file.path(datapath, 'ppiGraph.rda')  # undirected edge PPI
dppipath <- file.path(datapath, 'DppiGraph.rda')  # directed edge PPI
if(!file.exists(ppipath)) {
  print('ppiGraph does not exist, now creating ppiGraph start')
  cons_ppi(datapath, gdacpath)
}
load(file.path(ppipath))
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

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna", "meth")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                                method = "DRW", pranking = "t-test", classifier = "glm",
                                nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM, file=file.path('data/model/res_pa_GM.RData'))

summary(res_pa_GM)
print(res_pa_GM$results)
print(res_pa_GM$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GM, profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa(Pathway_Graph)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                                 method = "DRW", pranking = "t-test", classifier = "glm",
                                 nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR, file=file.path('data/model/res_pa_GMR.RData'))

summary(res_pa_GMR)
print(res_pa_GMR$results)
print(res_pa_GMR$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR, profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                                 method = "DRW", pranking = "t-test", classifier = "glm",
                                 nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP, file=file.path('data/model/res_pa_GMP.RData'))

summary(res_pa_GMP)
print(res_pa_GMP$results)
print(res_pa_GMP$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP, profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa(diffused_Pathway_Graph)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d <- fit.classification(y=y, samples = samples, respath = respath, profile_name = profile_name,
                                   method = "DRW", pranking = "t-test", classifier = "glm",
                                   nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d, file=file.path('data/model/res_pa_GMR_d.RData'))

summary(res_pa_GMR_d)
print(res_pa_GMR_d$results)
print(res_pa_GMR_d$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_d, profile_name=profile_name, method="DRW", respath=respath)

# plot
xlabs <- c("GM", "GMR", "GMR_d", "GMP")
res_models <- list(res_pa_GM, res_pa_GMR, res_pa_GMR_d, res_pa_GMP)

perf_boxplot(xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM$results$Accuracy)