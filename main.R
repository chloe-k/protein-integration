
########################### package load ############################
# source("https://bioconductor.org/biocLite.R")

# biocLite("KEGGgraph")
# biocLite("igraph")
# install.packages("ggplot2")
# biocLite("annotate")
# biocLite("org.Hs.eg.db")
# biocLite("diffusr")
# biocLite("DESeq2")
# biocLite("Matrix")
# biocLite("stringr")
# biocLite("caret")
# biocLite("e1071")
# biocLite("randomForest")
# biocLite("KEGG.db")
# biocLite("KEGGREST")

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


# ppipath <- file.path(datapath, 'ppiGraph.rda')  # undirected edge PPI
# dppipath <- file.path(datapath, 'DppiGraph.rda')  # directed edge PPI
#ppipath <- file.path(datapath, 'ppiGraph(Entrez).rda')  # undirected edge PPI
dppipath <- file.path(datapath, 'DppiGraph(Entrez).rda')  # directed edge PPI
# dppipath <- file.path(datapath, 'DppiGraph_rdc.rda')  # directed edge PPI
if(!file.exists(dppipath)) {
  print('ppiGraph does not exist, now creating ppiGraph start')
  cons_ppi(datapath, gdacpath, rppa)
}
#load(file.path(ppipath))
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

#######################################Result 19##########################################################
res_models <- list()
load('data/model/res_pa_GM_RF_13.RData')
load('data/model/res_pa_GMR_RF_13.RData')
res_models <- c(res_models, list(res_pa_GM_RF_13, res_pa_GMR_RF_13))

###########################################################  Prob = 0.01  #######################################################
prob <- 0.01
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.01 <- fit.classification(y=y, samples = samples, id = "result19", datapath = datapath, respath = respath, profile_name = profile_name,
                                           method = "DRW", pranking = "t-test", classifier = "rf",
                                           nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.01, file=file.path('data/model/res_pa_GMR_d_19_0.01.RData'))

summary(res_pa_GMR_d_19_0.01)
print(res_pa_GMR_d_19_0.01$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.01, id = "result19", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.01 <- fit.classification(y=y, samples = samples, id = "result19", datapath = datapath, respath = respath, profile_name = profile_name,
                                         method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.01, file=file.path('data/model/res_pa_GMP_19_0.01.RData'))

summary(res_pa_GMP_19_0.01)
print(res_pa_GMP_19_0.01$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.01, id = "result19", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.01, res_pa_GMP_19_0.01))

###########################################################  Prob = 0.05  #######################################################
prob <- 0.05
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.05", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.05 <- fit.classification(y=y, samples = samples, id = "result19_0.05", datapath = datapath, respath = respath, profile_name = profile_name,
                                           method = "DRW", pranking = "t-test", classifier = "rf",
                                           nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.05, file=file.path('data/model/res_pa_GMR_d_19_0.05.RData'))

summary(res_pa_GMR_d_19_0.05)
print(res_pa_GMR_d_19_0.05$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.05, id = "result19_0.05", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.05", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.05 <- fit.classification(y=y, samples = samples, id = "result19_0.05", datapath = datapath, respath = respath, profile_name = profile_name,
                                         method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.05, file=file.path('data/model/res_pa_GMP_19_0.05.RData'))

summary(res_pa_GMP_19_0.05)
print(res_pa_GMP_19_0.05$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.05, id = "result19_0.05", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.05, res_pa_GMP_19_0.05))

###########################################################  Prob = 0.1  #######################################################
prob <- 0.1
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.1", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.1 <- fit.classification(y=y, samples = samples, id = "result19_0.1", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.1, file=file.path('data/model/res_pa_GMR_d_19_0.1.RData'))

summary(res_pa_GMR_d_19_0.1)
print(res_pa_GMR_d_19_0.1$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.1, id = "result19_0.1", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.1", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.1 <- fit.classification(y=y, samples = samples, id = "result19_0.1", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.1, file=file.path('data/model/res_pa_GMP_19_0.1.RData'))

summary(res_pa_GMP_19_0.1)
print(res_pa_GMP_19_0.1$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.1, id = "result19_0.1", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.1, res_pa_GMP_19_0.1))

###########################################################  Prob = 0.2  #######################################################
prob <- 0.2
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.2", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.2 <- fit.classification(y=y, samples = samples, id = "result19_0.2", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.2, file=file.path('data/model/res_pa_GMR_d_19_0.2.RData'))

summary(res_pa_GMR_d_19_0.2)
print(res_pa_GMR_d_19_0.2$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.2, id = "result19_0.2", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.2", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.2 <- fit.classification(y=y, samples = samples, id = "result19_0.2", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.2, file=file.path('data/model/res_pa_GMP_19_0.2.RData'))

summary(res_pa_GMP_19_0.2)
print(res_pa_GMP_19_0.2$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.2, id = "result19_0.2", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.2, res_pa_GMP_19_0.2))

###########################################################  Prob = 0.3  #######################################################
prob <- 0.3
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.3", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.3 <- fit.classification(y=y, samples = samples, id = "result19_0.3", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.3, file=file.path('data/model/res_pa_GMR_d_19_0.3.RData'))

summary(res_pa_GMR_d_19_0.3)
print(res_pa_GMR_d_19_0.3$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.3, id = "result19_0.3", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.3", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.3 <- fit.classification(y=y, samples = samples, id = "result19_0.3", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.3, file=file.path('data/model/res_pa_GMP_19_0.3.RData'))

summary(res_pa_GMP_19_0.3)
print(res_pa_GMP_19_0.3$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.3, id = "result19_0.3", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.3, res_pa_GMP_19_0.3))

###########################################################  Prob = 0.4  #######################################################
prob <- 0.4
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.4", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.4 <- fit.classification(y=y, samples = samples, id = "result19_0.4", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.4, file=file.path('data/model/res_pa_GMR_d_19_0.4.RData'))

summary(res_pa_GMR_d_19_0.4)
print(res_pa_GMR_d_19_0.4$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.4, id = "result19_0.4", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.4", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.4 <- fit.classification(y=y, samples = samples, id = "result19_0.4", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.4, file=file.path('data/model/res_pa_GMP_19_0.4.RData'))

summary(res_pa_GMP_19_0.4)
print(res_pa_GMP_19_0.4$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.4, id = "result19_0.4", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.4, res_pa_GMP_19_0.4))

###########################################################  Prob = 0.5  #######################################################
prob <- 0.5
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.5", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.5 <- fit.classification(y=y, samples = samples, id = "result19_0.5", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.5, file=file.path('data/model/res_pa_GMR_d_19_0.5.RData'))

summary(res_pa_GMR_d_19_0.5)
print(res_pa_GMR_d_19_0.5$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.5, id = "result19_0.5", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.5", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.5 <- fit.classification(y=y, samples = samples, id = "result19_0.5", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.5, file=file.path('data/model/res_pa_GMP_19_0.5.RData'))

summary(res_pa_GMP_19_0.5)
print(res_pa_GMP_19_0.5$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.5, id = "result19_0.5", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.5, res_pa_GMP_19_0.5))

###########################################################  Prob = 0.6  #######################################################
prob <- 0.6
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.6", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.6 <- fit.classification(y=y, samples = samples, id = "result19_0.6", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.6, file=file.path('data/model/res_pa_GMR_d_19_0.6.RData'))

summary(res_pa_GMR_d_19_0.6)
print(res_pa_GMR_d_19_0.6$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.6, id = "result19_0.6", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.6", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.6 <- fit.classification(y=y, samples = samples, id = "result19_0.6", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.6, file=file.path('data/model/res_pa_GMP_19_0.6.RData'))

summary(res_pa_GMP_19_0.6)
print(res_pa_GMP_19_0.6$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.6, id = "result19_0.6", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.6, res_pa_GMP_19_0.6))

###########################################################  Prob = 0.7  #######################################################
prob <- 0.7
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.7", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.7 <- fit.classification(y=y, samples = samples, id = "result19_0.7", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.7, file=file.path('data/model/res_pa_GMR_d_19_0.7.RData'))

summary(res_pa_GMR_d_19_0.7)
print(res_pa_GMR_d_19_0.7$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.7, id = "result19_0.7", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.7", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.7 <- fit.classification(y=y, samples = samples, id = "result19_0.7", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.7, file=file.path('data/model/res_pa_GMP_19_0.7.RData'))

summary(res_pa_GMP_19_0.7)
print(res_pa_GMP_19_0.7$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.7, id = "result19_0.7", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.7, res_pa_GMP_19_0.7))


###########################################################  Prob = 0.8  #######################################################
prob <- 0.8
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.8", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.8 <- fit.classification(y=y, samples = samples, id = "result19_0.8", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.8, file=file.path('data/model/res_pa_GMR_d_19_0.8.RData'))

summary(res_pa_GMR_d_19_0.8)
print(res_pa_GMR_d_19_0.8$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.8, id = "result19_0.8", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.8", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.8 <- fit.classification(y=y, samples = samples, id = "result19_0.8", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.8, file=file.path('data/model/res_pa_GMP_19_0.8.RData'))

summary(res_pa_GMP_19_0.8)
print(res_pa_GMP_19_0.8$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.8, id = "result19_0.8", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.8, res_pa_GMP_19_0.8))

###########################################################  Prob = 0.9  #######################################################
prob <- 0.9
#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.9", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_19_0.9 <- fit.classification(y=y, samples = samples, id = "result19_0.9", datapath = datapath, respath = respath, profile_name = profile_name,
                                          method = "DRW", pranking = "t-test", classifier = "rf",
                                          nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_19_0.9, file=file.path('data/model/res_pa_GMR_d_19_0.9.RData'))

summary(res_pa_GMR_d_19_0.9)
print(res_pa_GMR_d_19_0.9$results)

write.SigFeatures(res_fit=res_pa_GMR_d_19_0.9, id = "result19_0.9", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result19_0.9", prob = prob,
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_19_0.9 <- fit.classification(y=y, samples = samples, id = "result19_0.9", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_19_0.9, file=file.path('data/model/res_pa_GMP_19_0.9.RData'))

summary(res_pa_GMP_19_0.9)
print(res_pa_GMP_19_0.9$results)

write.SigFeatures(res_fit=res_pa_GMP_19_0.9, id = "result19_0.9", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_19_0.9, res_pa_GMP_19_0.9))


# plot
title <- c("Result 19")
xlabs <- c("GM", "GMR", "GMR_d_0.01", "GMP_0.01", "GMR_d_0.05", "GMP_0.05", "GMR_d_0.1", "GMP_0.1", "GMR_d_0.2", "GMP_0.2", "GMR_d_0.3", "GMP_0.3", "GMR_d_0.4", "GMP_0.4", "GMR_d_0.5", "GMP_0.5", "GMR_d_0.6", "GMP_0.6", "GMR_d_0.7", "GMP_0.7", "GMR_d_0.8", "GMP_0.8", "GMR_d_0.9", "GMP_0.9")
#res_models <- list(res_pa_GM_13, res_pa_GMR_13, res_pa_GMR_d_19_0.01, res_pa_GMP_19_0.01, res_pa_GMR_d_19_0.05, res_pa_GMP_19_0.05, res_pa_GMR_d_19_0.1, res_pa_GMP_19_0.1)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_RF_13$results$Accuracy[1])
