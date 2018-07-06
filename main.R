
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

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result9",
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_RF_9 <- fit.classification(y=y, samples = samples, id = "result9", datapath = datapath, respath = respath, profile_name = profile_name,
                                     method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_RF_9, file=file.path('data/model/res_pa_GM_RF_9.RData'))

summary(res_pa_GM_RF_9)
print(res_pa_GM_RF_9$results)

write.SigFeatures(res_fit=res_pa_GM_RF_9, id = "result9", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result9",
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_RF_9 <- fit.classification(y=y, samples = samples, id = "result9", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_RF_9, file=file.path('data/model/res_pa_GMR_RF_9.RData'))

summary(res_pa_GMR_RF_9)
print(res_pa_GMR_RF_9$results)

write.SigFeatures(res_fit=res_pa_GMR_RF_9, id = "result9", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result9",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_RF_9 <- fit.classification(y=y, samples = samples, id = "result9", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_RF_9, file=file.path('data/model/res_pa_GMR_d_RF_9.RData'))

summary(res_pa_GMR_d_RF_9)
print(res_pa_GMR_d_RF_9$results)

write.SigFeatures(res_fit=res_pa_GMR_d_RF_9, id = "result9", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result9",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_RF_9 <- fit.classification(y=y, samples = samples, id = "result9", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_RF_9, file=file.path('data/model/res_pa_GMP_RF_9.RData'))

summary(res_pa_GMP_RF_9)
print(res_pa_GMP_RF_9$results)

write.SigFeatures(res_fit=res_pa_GMP_RF_9, id = "result9", profile_name=profile_name, method="DRW", respath=respath)


# plot
title <- c("Result 9")
xlabs <- c("GM", "GMR", "GMR_d", "GMP")
res_models <- list(res_pa_GM_RF_9, res_pa_GMR_RF_9, res_pa_GMR_d_RF_9, res_pa_GMP_RF_9)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_RF_9$results$Accuracy[1])


###################################Result 16#####################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result16",
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_RF_16 <- fit.classification(y=y, samples = samples, id = "result16", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_RF_16, file=file.path('data/model/res_pa_GM_RF_16.RData'))

summary(res_pa_GM_RF_16)
print(res_pa_GM_RF_16$results)
print(res_pa_GM_RF_16$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GM_RF_16, id = "result16", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result16",
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_RF_16 <- fit.classification(y=y, samples = samples, id = "result16", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_RF_16, file=file.path('data/model/res_pa_GMR_RF_16.RData'))

summary(res_pa_GMR_RF_16)
print(res_pa_GMR_RF_16$results)
print(res_pa_GMR_RF_16$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_RF_16, id = "result16", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result16",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_RF_16 <- fit.classification(y=y, samples = samples, id = "result16", datapath = datapath, respath = respath, profile_name = profile_name,
                                         method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_RF_16, file=file.path('data/model/res_pa_GMR_d_RF_16.RData'))

summary(res_pa_GMR_d_RF_16)
print(res_pa_GMR_d_RF_16$results)
print(res_pa_GMR_d_RF_16$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_d_RF_16, id = "result16", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id = "result16",
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_RF_16 <- fit.classification(y=y, samples = samples, id = "result16", datapath = datapath, respath = respath, profile_name = profile_name,
                                       method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_RF_16, file=file.path('data/model/res_pa_GMP_RF_16.RData'))

summary(res_pa_GMP_RF_16)
print(res_pa_GMP_RF_16$results)
print(res_pa_GMP_RF_16$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_RF_16, id = "result16", profile_name=profile_name, method="DRW", respath=respath)


# plot
title <- c("Result 16")
xlabs <- c("GM", "GMR", "GMR_d", "GMP")
res_models <- list(res_pa_GM_RF_16, res_pa_GMR_RF_16, res_pa_GMR_d_RF_16, res_pa_GMP_RF_16)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_RF_16$results$Accuracy[1])


###################################Result 10_all#####################################################

#---------------------------------Result 10_1---------------------------------#
# PPI relation : interacts-with
dppipath <- file.path(datapath, 'DppiGraph_1.rda')
load(file.path(dppipath))

p <- DppiGraph_1
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_1", 
               pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_10_1 <- fit.classification(y=y, samples = samples, id="result10_1", datapath = datapath, respath = respath, profile_name = profile_name,
                                     method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_10_1, file=file.path('data/model/res_pa_GM_10_1.RData'))

summary(res_pa_GM_10_1)
print(res_pa_GM_10_1$results)
print(res_pa_GM_10_1$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GM_10_1, id="result10_1", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_1", 
               pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_10_1 <- fit.classification(y=y, samples = samples, id="result10_1", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_10_1, file=file.path('data/model/res_pa_GMR_10_1.RData'))

summary(res_pa_GMR_10_1)
print(res_pa_GMR_10_1$results)
print(res_pa_GMR_10_1$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_10_1, id="result10_1", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_1", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_10_1 <- fit.classification(y=y, samples = samples, id="result10_1", datapath = datapath, respath = respath, profile_name = profile_name,
                                        method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_10_1, file=file.path('data/model/res_pa_GMR_d_10_1.RData'))

summary(res_pa_GMR_d_10_1)
print(res_pa_GMR_d_10_1$results)
print(res_pa_GMR_d_10_1$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMR_d_10_1, id="result10_1", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_1", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_1 <- fit.classification(y=y, samples = samples, id="result10_1", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_1, file=file.path('data/model/res_pa_GMP_10_1.RData'))

summary(res_pa_GMP_10_1)
print(res_pa_GMP_10_1$results)
print(res_pa_GMP_10_1$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_1, id="result10_1", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_2---------------------------------#
# PPI relation : controls-phosphorylation-of
dppipath <- file.path(datapath, 'DppiGraph_2.rda')
load(file.path(dppipath))

p <- DppiGraph_2
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_2", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_2 <- fit.classification(y=y, samples = samples, id="result10_2", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_2, file=file.path('data/model/res_pa_GMP_10_2.RData'))

summary(res_pa_GMP_10_2)
print(res_pa_GMP_10_2$results)
print(res_pa_GMP_10_2$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_2, id="result10_2", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_3---------------------------------#
# PPI relation : catalysis-precedes
dppipath <- file.path(datapath, 'DppiGraph_3.rda')
load(file.path(dppipath))

p <- DppiGraph_3
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_3", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_3 <- fit.classification(y=y, samples = samples, id="result10_3", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_3, file=file.path('data/model/res_pa_GMP_10_3.RData'))

summary(res_pa_GMP_10_3)
print(res_pa_GMP_10_3$results)
print(res_pa_GMP_10_3$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_3, id="result10_3", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_4---------------------------------#
# PPI relation : controls-expression-of
dppipath <- file.path(datapath, 'DppiGraph_4.rda')
load(file.path(dppipath))

p <- DppiGraph_4
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_4", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_4 <- fit.classification(y=y, samples = samples, id="result10_4", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_4, file=file.path('data/model/res_pa_GMP_10_4.RData'))

summary(res_pa_GMP_10_4)
print(res_pa_GMP_10_4$results)
print(res_pa_GMP_10_4$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_4, id="result10_4", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_5---------------------------------#
# PPI relation : controls-state-change-of
dppipath <- file.path(datapath, 'DppiGraph_5.rda')
load(file.path(dppipath))

p <- DppiGraph_5
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_5", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_5 <- fit.classification(y=y, samples = samples, id="result10_5", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_5, file=file.path('data/model/res_pa_GMP_10_5.RData'))

summary(res_pa_GMP_10_5)
print(res_pa_GMP_10_5$results)
print(res_pa_GMP_10_5$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_5, id="result10_5", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_6---------------------------------#
# PPI relation : in-complex-with
dppipath <- file.path(datapath, 'DppiGraph_6.rda')
load(file.path(dppipath))

p <- DppiGraph_6
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_6", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_6 <- fit.classification(y=y, samples = samples, id="result10_6", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_6, file=file.path('data/model/res_pa_GMP_10_6.RData'))

summary(res_pa_GMP_10_6)
print(res_pa_GMP_10_6$results)
print(res_pa_GMP_10_6$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_6, id="result10_6", profile_name=profile_name, method="DRW", respath=respath)

#---------------------------------Result 10_7---------------------------------#
# PPI relation : controls-transport-of
dppipath <- file.path(datapath, 'DppiGraph_7.rda')
load(file.path(dppipath))

p <- DppiGraph_7
V(p)$name <-paste("p",V(p)$name,sep="")

#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, id="result10_7", 
               pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_10_7 <- fit.classification(y=y, samples = samples, id="result10_7", datapath = datapath, respath = respath, profile_name = profile_name,
                                      method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_10_7, file=file.path('data/model/res_pa_GMP_10_7.RData'))

summary(res_pa_GMP_10_7)
print(res_pa_GMP_10_7$results)
print(res_pa_GMP_10_7$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP_10_7, id="result10_7", profile_name=profile_name, method="DRW", respath=respath)



# plot
title <- c("Result 10_all")
xlabs <- c("GM", "GMR", "GMR_d", "GMP_1", "GMP_2", "GMP_3", "GMP_4", "GMP_5", "GMP_6", "GMP_7")
res_models <- list(res_pa_GM_10_1, res_pa_GMR_10_1, res_pa_GMR_d_10_1, res_pa_GMP_10_1, res_pa_GMP_10_2, res_pa_GMP_10_3, res_pa_GMP_10_4, res_pa_GMP_10_5, res_pa_GMP_10_6, res_pa_GMP_10_7)

perf_boxplot(title, xlabs, res_models, perf_min = 0.5, perf_max = 0.9, res_pa_GM_10_1$results$Accuracy[1])