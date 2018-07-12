
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




################################## Result 18_all ############################################################

res_models <- list()
res_gm <- list()
res_gmr <- list()
res_gmr_d <- list()
res_gmp <- list()


################################################### Result18_1: prob = 0.2, Gamma = 0.2  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_1_GM", prob = 0.2, Gamma = 0.2, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_1 <- fit.classification(y=y, samples = samples, id = "result18_1_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_1, file=file.path('data/model/res_pa_GM_18_1.RData'))

summary(res_pa_GM_18_1)

write.SigFeatures(res_fit=res_pa_GM_18_1, id = "result18_1_GM", profile_name=profile_name, method="DRW", respath=respath)

res_gm <- c(res_gm, list(res_pa_GM_18_1))

#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_1_GMR", prob = 0.2, Gamma = 0.2, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_1 <- fit.classification(y=y, samples = samples, id = "result18_1_GMR", datapath = datapath, respath = respath, 
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_1, file=file.path('data/model/res_pa_GMR_18_1.RData'))

summary(res_pa_GMR_18_1)

write.SigFeatures(res_fit=res_pa_GMR_18_1, id = "result18_1_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_1_GMR_d", prob = 0.2, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_1 <- fit.classification(y=y, samples = samples, id = "result18_1_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_1, file=file.path('data/model/res_pa_GMR_d_18_1.RData'))

summary(res_pa_GMR_d_18_1)

write.SigFeatures(res_fit=res_pa_GMR_d_18_1, id = "result18_1_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_1_GMP", prob = 0.2, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_1 <- fit.classification(y=y, samples = samples, id = "result18_1_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_1, file=file.path('data/model/res_pa_GMP_18_1.RData'))

summary(res_pa_GMP_18_1)

write.SigFeatures(res_fit=res_pa_GMP_18_1, id = "result18_1_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_1, res_pa_GMR_18_1, res_pa_GMR_d_18_1, res_pa_GMP_18_1))


################################################### Result18_2: prob = 0.2, Gamma = 0.4  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_2_GM", prob = 0.2, Gamma = 0.4, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_2 <- fit.classification(y=y, samples = samples, id = "result18_2_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_2, file=file.path('data/model/res_pa_GM_18_2.RData'))

summary(res_pa_GM_18_2)

write.SigFeatures(res_fit=res_pa_GM_18_2, id = "result18_2_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_2_GMR", prob = 0.2, Gamma = 0.4, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_2 <- fit.classification(y=y, samples = samples, id = "result18_2_GMR", datapath = datapath, respath = respath, 
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_2, file=file.path('data/model/res_pa_GMR_18_2.RData'))

summary(res_pa_GMR_18_2)

write.SigFeatures(res_fit=res_pa_GMR_18_2, id = "result18_2_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_2_GMR_d", prob = 0.2, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_2 <- fit.classification(y=y, samples = samples, id = "result18_2_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_2, file=file.path('data/model/res_pa_GMR_d_18_2.RData'))

summary(res_pa_GMR_d_18_2)

write.SigFeatures(res_fit=res_pa_GMR_d_18_2, id = "result18_2_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_2_GMP", prob = 0.2, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_2 <- fit.classification(y=y, samples = samples, id = "result18_2_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_2, file=file.path('data/model/res_pa_GMP_18_2.RData'))

summary(res_pa_GMP_18_2)

write.SigFeatures(res_fit=res_pa_GMP_18_2, id = "result18_2_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_2, res_pa_GMR_18_2, res_pa_GMR_d_18_2, res_pa_GMP_18_2))


################################################### Result18_3: prob = 0.2, Gamma = 0.6  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_3_GM", prob = 0.2, Gamma = 0.6, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_3 <- fit.classification(y=y, samples = samples, id = "result18_3_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_3, file=file.path('data/model/res_pa_GM_18_3.RData'))

summary(res_pa_GM_18_3)

write.SigFeatures(res_fit=res_pa_GM_18_3, id = "result18_3_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_3_GMR", prob = 0.2, Gamma = 0.6, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_3 <- fit.classification(y=y, samples = samples, id = "result18_3_GMR", datapath = datapath, respath = respath, 
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_3, file=file.path('data/model/res_pa_GMR_18_3.RData'))

summary(res_pa_GMR_18_3)

write.SigFeatures(res_fit=res_pa_GMR_18_3, id = "result18_3_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_3_GMR_d", prob = 0.2, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_3 <- fit.classification(y=y, samples = samples, id = "result18_3_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_3, file=file.path('data/model/res_pa_GMR_d_18_3.RData'))

summary(res_pa_GMR_d_18_3)

write.SigFeatures(res_fit=res_pa_GMR_d_18_3, id = "result18_3_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_3_GMP", prob = 0.2, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_3 <- fit.classification(y=y, samples = samples, id = "result18_3_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_3, file=file.path('data/model/res_pa_GMP_18_3.RData'))

summary(res_pa_GMP_18_3)

write.SigFeatures(res_fit=res_pa_GMP_18_3, id = "result18_3_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_3, res_pa_GMR_18_3, res_pa_GMR_d_18_3, res_pa_GMP_18_3))

################################################### Result18_4: prob = 0.2, Gamma = 0.8  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_4_GM", prob = 0.2, Gamma = 0.8, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_4 <- fit.classification(y=y, samples = samples, id = "result18_4_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_4, file=file.path('data/model/res_pa_GM_18_4.RData'))

summary(res_pa_GM_18_4)

write.SigFeatures(res_fit=res_pa_GM_18_4, id = "result18_4_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_4_GMR", prob = 0.2, Gamma = 0.8, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_4 <- fit.classification(y=y, samples = samples, id = "result18_4_GMR", datapath = datapath, respath = respath, 
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_4, file=file.path('data/model/res_pa_GMR_18_4.RData'))

summary(res_pa_GMR_18_4)

write.SigFeatures(res_fit=res_pa_GMR_18_4, id = "result18_4_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_4_GMR_d", prob = 0.2, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_4 <- fit.classification(y=y, samples = samples, id = "result18_4_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_4, file=file.path('data/model/res_pa_GMR_d_18_4.RData'))

summary(res_pa_GMR_d_18_4)

write.SigFeatures(res_fit=res_pa_GMR_d_18_4, id = "result18_4_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_4_GMP", prob = 0.2, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_4 <- fit.classification(y=y, samples = samples, id = "result18_4_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_4, file=file.path('data/model/res_pa_GMP_18_4.RData'))

summary(res_pa_GMP_18_4)

write.SigFeatures(res_fit=res_pa_GMP_18_4, id = "result18_4_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_4, res_pa_GMR_18_4, res_pa_GMR_d_18_4, res_pa_GMP_18_4))


################################################### Result18_5: prob = 0.4, Gamma = 0.2  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_5_GM", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_5 <- fit.classification(y=y, samples = samples, id = "result18_5_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_5, file=file.path('data/model/res_pa_GM_18_5.RData'))

summary(res_pa_GM_18_5)

write.SigFeatures(res_fit=res_pa_GM_18_5, id = "result18_5_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_5_GMR", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_5 <- fit.classification(y=y, samples = samples, id = "result18_5_GMR", datapath = datapath, respath = respath, 
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_5, file=file.path('data/model/res_pa_GMR_18_5.RData'))

summary(res_pa_GMR_18_5)

write.SigFeatures(res_fit=res_pa_GMR_18_5, id = "result18_5_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_5_GMR_d", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_5 <- fit.classification(y=y, samples = samples, id = "result18_5_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_5, file=file.path('data/model/res_pa_GMR_d_18_5.RData'))

summary(res_pa_GMR_d_18_5)

write.SigFeatures(res_fit=res_pa_GMR_d_18_5, id = "result18_5_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_5_GMP", prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_5 <- fit.classification(y=y, samples = samples, id = "result18_5_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_5, file=file.path('data/model/res_pa_GMP_18_5.RData'))

summary(res_pa_GMP_18_5)

write.SigFeatures(res_fit=res_pa_GMP_18_5, id = "result18_5_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_5, res_pa_GMR_18_5, res_pa_GMR_d_18_5, res_pa_GMP_18_5))


################################################### Result18_6: prob = 0.4, Gamma = 0.4  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_6_GM", prob = 0.4, Gamma = 0.4, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_6 <- fit.classification(y=y, samples = samples, id = "result18_6_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_6, file=file.path('data/model/res_pa_GM_18_6.RData'))

summary(res_pa_GM_18_6)

write.SigFeatures(res_fit=res_pa_GM_18_6, id = "result18_6_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_6_GMR", prob = 0.4, Gamma = 0.4, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_6 <- fit.classification(y=y, samples = samples, id = "result18_6_GMR", datapath = datapath, respath = respath, 
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_6, file=file.path('data/model/res_pa_GMR_18_6.RData'))

summary(res_pa_GMR_18_6)

write.SigFeatures(res_fit=res_pa_GMR_18_6, id = "result18_6_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_6_GMR_d", prob = 0.4, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_6 <- fit.classification(y=y, samples = samples, id = "result18_6_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_6, file=file.path('data/model/res_pa_GMR_d_18_6.RData'))

summary(res_pa_GMR_d_18_6)

write.SigFeatures(res_fit=res_pa_GMR_d_18_6, id = "result18_6_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_6_GMP", prob = 0.4, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_6 <- fit.classification(y=y, samples = samples, id = "result18_6_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_6, file=file.path('data/model/res_pa_GMP_18_6.RData'))

summary(res_pa_GMP_18_6)

write.SigFeatures(res_fit=res_pa_GMP_18_6, id = "result18_6_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_6, res_pa_GMR_18_6, res_pa_GMR_d_18_6, res_pa_GMP_18_6))


################################################### Result18_7: prob = 0.4, Gamma = 0.6  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_7_GM", prob = 0.4, Gamma = 0.6, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_7 <- fit.classification(y=y, samples = samples, id = "result18_7_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_7, file=file.path('data/model/res_pa_GM_18_7.RData'))

summary(res_pa_GM_18_7)

write.SigFeatures(res_fit=res_pa_GM_18_7, id = "result18_7_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_7_GMR", prob = 0.4, Gamma = 0.6, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_7 <- fit.classification(y=y, samples = samples, id = "result18_7_GMR", datapath = datapath, respath = respath, 
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_7, file=file.path('data/model/res_pa_GMR_18_7.RData'))

summary(res_pa_GMR_18_7)

write.SigFeatures(res_fit=res_pa_GMR_18_7, id = "result18_7_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_7_GMR_d", prob = 0.4, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_7 <- fit.classification(y=y, samples = samples, id = "result18_7_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_7, file=file.path('data/model/res_pa_GMR_d_18_7.RData'))

summary(res_pa_GMR_d_18_7)

write.SigFeatures(res_fit=res_pa_GMR_d_18_7, id = "result18_7_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_7_GMP", prob = 0.4, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_7 <- fit.classification(y=y, samples = samples, id = "result18_7_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_7, file=file.path('data/model/res_pa_GMP_18_7.RData'))

summary(res_pa_GMP_18_7)

write.SigFeatures(res_fit=res_pa_GMP_18_7, id = "result18_7_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_7, res_pa_GMR_18_7, res_pa_GMR_d_18_7, res_pa_GMP_18_7))


################################################### Result18_8: prob = 0.4, Gamma = 0.8  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_8_GM", prob = 0.4, Gamma = 0.8, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_8 <- fit.classification(y=y, samples = samples, id = "result18_8_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_8, file=file.path('data/model/res_pa_GM_18_8.RData'))

summary(res_pa_GM_18_8)

write.SigFeatures(res_fit=res_pa_GM_18_8, id = "result18_8_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_8_GMR", prob = 0.4, Gamma = 0.8, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_8 <- fit.classification(y=y, samples = samples, id = "result18_8_GMR", datapath = datapath, respath = respath, 
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_8, file=file.path('data/model/res_pa_GMR_18_8.RData'))

summary(res_pa_GMR_18_8)

write.SigFeatures(res_fit=res_pa_GMR_18_8, id = "result18_8_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_8_GMR_d", prob = 0.4, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_8 <- fit.classification(y=y, samples = samples, id = "result18_8_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_8, file=file.path('data/model/res_pa_GMR_d_18_8.RData'))

summary(res_pa_GMR_d_18_8)

write.SigFeatures(res_fit=res_pa_GMR_d_18_8, id = "result18_8_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_8_GMP", prob = 0.4, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_8 <- fit.classification(y=y, samples = samples, id = "result18_8_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_8, file=file.path('data/model/res_pa_GMP_18_8.RData'))

summary(res_pa_GMP_18_8)

write.SigFeatures(res_fit=res_pa_GMP_18_8, id = "result18_8_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_8, res_pa_GMR_18_8, res_pa_GMR_d_18_8, res_pa_GMP_18_8))


################################################### Result18_9: prob = 0.6, Gamma = 0.2  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_9_GM", prob = 0.6, Gamma = 0.2, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_9 <- fit.classification(y=y, samples = samples, id = "result18_9_GM", datapath = datapath, respath = respath,
                                     profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                     nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_9, file=file.path('data/model/res_pa_GM_18_9.RData'))

summary(res_pa_GM_18_9)

write.SigFeatures(res_fit=res_pa_GM_18_9, id = "result18_9_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_9_GMR", prob = 0.6, Gamma = 0.2, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_9 <- fit.classification(y=y, samples = samples, id = "result18_9_GMR", datapath = datapath, respath = respath, 
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_9, file=file.path('data/model/res_pa_GMR_18_9.RData'))

summary(res_pa_GMR_18_9)

write.SigFeatures(res_fit=res_pa_GMR_18_9, id = "result18_9_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_9_GMR_d", prob = 0.6, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_9 <- fit.classification(y=y, samples = samples, id = "result18_9_GMR_d", datapath = datapath, respath = respath,
                                        profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                        nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_9, file=file.path('data/model/res_pa_GMR_d_18_9.RData'))

summary(res_pa_GMR_d_18_9)

write.SigFeatures(res_fit=res_pa_GMR_d_18_9, id = "result18_9_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_9_GMP", prob = 0.6, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_9 <- fit.classification(y=y, samples = samples, id = "result18_9_GMP", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_9, file=file.path('data/model/res_pa_GMP_18_9.RData'))

summary(res_pa_GMP_18_9)

write.SigFeatures(res_fit=res_pa_GMP_18_9, id = "result18_9_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_9, res_pa_GMR_18_9, res_pa_GMR_d_18_9, res_pa_GMP_18_9))


################################################### Result18_10: prob = 0.6, Gamma = 0.4  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_10_GM", prob = 0.6, Gamma = 0.4, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_10 <- fit.classification(y=y, samples = samples, id = "result18_10_GM", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_10, file=file.path('data/model/res_pa_GM_18_10.RData'))

summary(res_pa_GM_18_10)

write.SigFeatures(res_fit=res_pa_GM_18_10, id = "result18_10_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_10_GMR", prob = 0.6, Gamma = 0.4, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_10 <- fit.classification(y=y, samples = samples, id = "result18_10_GMR", datapath = datapath, respath = respath, 
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_10, file=file.path('data/model/res_pa_GMR_18_10.RData'))

summary(res_pa_GMR_18_10)

write.SigFeatures(res_fit=res_pa_GMR_18_10, id = "result18_10_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_10_GMR_d", prob = 0.6, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_10 <- fit.classification(y=y, samples = samples, id = "result18_10_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_10, file=file.path('data/model/res_pa_GMR_d_18_10.RData'))

summary(res_pa_GMR_d_18_10)

write.SigFeatures(res_fit=res_pa_GMR_d_18_10, id = "result18_10_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_10_GMP", prob = 0.6, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_10 <- fit.classification(y=y, samples = samples, id = "result18_10_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_10, file=file.path('data/model/res_pa_GMP_18_10.RData'))

summary(res_pa_GMP_18_10)

write.SigFeatures(res_fit=res_pa_GMP_18_10, id = "result18_10_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_10, res_pa_GMR_18_10, res_pa_GMR_d_18_10, res_pa_GMP_18_10))


################################################### Result18_11: prob = 0.6, Gamma = 0.6  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_11_GM", prob = 0.6, Gamma = 0.6, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_11 <- fit.classification(y=y, samples = samples, id = "result18_11_GM", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_11, file=file.path('data/model/res_pa_GM_18_11.RData'))

summary(res_pa_GM_18_11)

write.SigFeatures(res_fit=res_pa_GM_18_11, id = "result18_11_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_11_GMR", prob = 0.6, Gamma = 0.6, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_11 <- fit.classification(y=y, samples = samples, id = "result18_11_GMR", datapath = datapath, respath = respath, 
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_11, file=file.path('data/model/res_pa_GMR_18_11.RData'))

summary(res_pa_GMR_18_11)

write.SigFeatures(res_fit=res_pa_GMR_18_11, id = "result18_11_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_11_GMR_d", prob = 0.6, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_11 <- fit.classification(y=y, samples = samples, id = "result18_11_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_11, file=file.path('data/model/res_pa_GMR_d_18_11.RData'))

summary(res_pa_GMR_d_18_11)

write.SigFeatures(res_fit=res_pa_GMR_d_18_11, id = "result18_11_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_11_GMP", prob = 0.6, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_11 <- fit.classification(y=y, samples = samples, id = "result18_11_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_11, file=file.path('data/model/res_pa_GMP_18_11.RData'))

summary(res_pa_GMP_18_11)

write.SigFeatures(res_fit=res_pa_GMP_18_11, id = "result18_11_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_11, res_pa_GMR_18_11, res_pa_GMR_d_18_11, res_pa_GMP_18_11))


################################################### Result18_12: prob = 0.6, Gamma = 0.8  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_12_GM", prob = 0.6, Gamma = 0.8, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_12 <- fit.classification(y=y, samples = samples, id = "result18_12_GM", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_12, file=file.path('data/model/res_pa_GM_18_12.RData'))

summary(res_pa_GM_18_12)

write.SigFeatures(res_fit=res_pa_GM_18_12, id = "result18_12_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_12_GMR", prob = 0.6, Gamma = 0.8, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_12 <- fit.classification(y=y, samples = samples, id = "result18_12_GMR", datapath = datapath, respath = respath, 
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_12, file=file.path('data/model/res_pa_GMR_18_12.RData'))

summary(res_pa_GMR_18_12)

write.SigFeatures(res_fit=res_pa_GMR_18_12, id = "result18_12_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_12_GMR_d", prob = 0.6, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_12 <- fit.classification(y=y, samples = samples, id = "result18_12_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_12, file=file.path('data/model/res_pa_GMR_d_18_12.RData'))

summary(res_pa_GMR_d_18_12)

write.SigFeatures(res_fit=res_pa_GMR_d_18_12, id = "result18_12_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_12_GMP", prob = 0.6, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_12 <- fit.classification(y=y, samples = samples, id = "result18_12_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_12, file=file.path('data/model/res_pa_GMP_18_12.RData'))

summary(res_pa_GMP_18_12)

write.SigFeatures(res_fit=res_pa_GMP_18_12, id = "result18_12_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_12, res_pa_GMR_18_12, res_pa_GMR_d_18_12, res_pa_GMP_18_12))


################################################### Result18_13: prob = 0.8, Gamma = 0.2  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_13_GM", prob = 0.8, Gamma = 0.2, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_13 <- fit.classification(y=y, samples = samples, id = "result18_13_GM", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_13, file=file.path('data/model/res_pa_GM_18_13.RData'))

summary(res_pa_GM_18_13)

write.SigFeatures(res_fit=res_pa_GM_18_13, id = "result18_13_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_13_GMR", prob = 0.8, Gamma = 0.2, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_13 <- fit.classification(y=y, samples = samples, id = "result18_13_GMR", datapath = datapath, respath = respath, 
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_13, file=file.path('data/model/res_pa_GMR_18_13.RData'))

summary(res_pa_GMR_18_13)

write.SigFeatures(res_fit=res_pa_GMR_18_13, id = "result18_13_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_13_GMR_d", prob = 0.8, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_13 <- fit.classification(y=y, samples = samples, id = "result18_13_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_13, file=file.path('data/model/res_pa_GMR_d_18_13.RData'))

summary(res_pa_GMR_d_18_13)

write.SigFeatures(res_fit=res_pa_GMR_d_18_13, id = "result18_13_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_13_GMP", prob = 0.8, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_13 <- fit.classification(y=y, samples = samples, id = "result18_13_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_13, file=file.path('data/model/res_pa_GMP_18_13.RData'))

summary(res_pa_GMP_18_13)

write.SigFeatures(res_fit=res_pa_GMP_18_13, id = "result18_13_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_13, res_pa_GMR_18_13, res_pa_GMR_d_18_13, res_pa_GMP_18_13))


################################################### Result18_14: prob = 0.8, Gamma = 0.4  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_14_GM", prob = 0.8, Gamma = 0.4, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_14 <- fit.classification(y=y, samples = samples, id = "result18_14_GM", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_14, file=file.path('data/model/res_pa_GM_18_14.RData'))

summary(res_pa_GM_18_14)

write.SigFeatures(res_fit=res_pa_GM_18_14, id = "result18_14_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_14_GMR", prob = 0.8, Gamma = 0.4, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_14 <- fit.classification(y=y, samples = samples, id = "result18_14_GMR", datapath = datapath, respath = respath, 
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_14, file=file.path('data/model/res_pa_GMR_18_14.RData'))

summary(res_pa_GMR_18_14)

write.SigFeatures(res_fit=res_pa_GMR_18_14, id = "result18_14_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_14_GMR_d", prob = 0.8, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_14 <- fit.classification(y=y, samples = samples, id = "result18_14_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_14, file=file.path('data/model/res_pa_GMR_d_18_14.RData'))

summary(res_pa_GMR_d_18_14)

write.SigFeatures(res_fit=res_pa_GMR_d_18_14, id = "result18_14_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_14_GMP", prob = 0.8, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_14 <- fit.classification(y=y, samples = samples, id = "result18_14_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_14, file=file.path('data/model/res_pa_GMP_18_14.RData'))

summary(res_pa_GMP_18_14)

write.SigFeatures(res_fit=res_pa_GMP_18_14, id = "result18_14_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_14, res_pa_GMR_18_14, res_pa_GMR_d_18_14, res_pa_GMP_18_14))


################################################### Result18_15: prob = 0.8, Gamma = 0.6  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_15_GM", prob = 0.8, Gamma = 0.6, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_15 <- fit.classification(y=y, samples = samples, id = "result18_15_GM", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_15, file=file.path('data/model/res_pa_GM_18_15.RData'))

summary(res_pa_GM_18_15)

write.SigFeatures(res_fit=res_pa_GM_18_15, id = "result18_15_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_15_GMR", prob = 0.8, Gamma = 0.6, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_15 <- fit.classification(y=y, samples = samples, id = "result18_15_GMR", datapath = datapath, respath = respath, 
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_15, file=file.path('data/model/res_pa_GMR_18_15.RData'))

summary(res_pa_GMR_18_15)

write.SigFeatures(res_fit=res_pa_GMR_18_15, id = "result18_15_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_15_GMR_d", prob = 0.8, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_15 <- fit.classification(y=y, samples = samples, id = "result18_15_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_15, file=file.path('data/model/res_pa_GMR_d_18_15.RData'))

summary(res_pa_GMR_d_18_15)

write.SigFeatures(res_fit=res_pa_GMR_d_18_15, id = "result18_15_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_15_GMP", prob = 0.8, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_15 <- fit.classification(y=y, samples = samples, id = "result18_15_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_15, file=file.path('data/model/res_pa_GMP_18_15.RData'))

summary(res_pa_GMP_18_15)

write.SigFeatures(res_fit=res_pa_GMP_18_15, id = "result18_15_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_15, res_pa_GMR_18_15, res_pa_GMR_d_18_15, res_pa_GMP_18_15))


################################################### Result18_16: prob = 0.8, Gamma = 0.8  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_16_GM", prob = 0.8, Gamma = 0.8, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_16 <- fit.classification(y=y, samples = samples, id = "result18_16_GM", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_16, file=file.path('data/model/res_pa_GM_18_16.RData'))

summary(res_pa_GM_18_16)

write.SigFeatures(res_fit=res_pa_GM_18_16, id = "result18_16_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_16_GMR", prob = 0.8, Gamma = 0.8, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_16 <- fit.classification(y=y, samples = samples, id = "result18_16_GMR", datapath = datapath, respath = respath, 
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_16, file=file.path('data/model/res_pa_GMR_18_16.RData'))

summary(res_pa_GMR_18_16)

write.SigFeatures(res_fit=res_pa_GMR_18_16, id = "result18_16_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_16_GMR_d", prob = 0.8, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_16 <- fit.classification(y=y, samples = samples, id = "result18_16_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_16, file=file.path('data/model/res_pa_GMR_d_18_16.RData'))

summary(res_pa_GMR_d_18_16)

write.SigFeatures(res_fit=res_pa_GMR_d_18_16, id = "result18_16_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_16_GMP", prob = 0.8, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_16 <- fit.classification(y=y, samples = samples, id = "result18_16_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_16, file=file.path('data/model/res_pa_GMP_18_16.RData'))

summary(res_pa_GMP_18_16)

write.SigFeatures(res_fit=res_pa_GMP_18_16, id = "result18_16_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_16, res_pa_GMR_18_16, res_pa_GMR_d_18_16, res_pa_GMP_18_16))


################################################### Result18_17: prob = 0, Gamma = 0  #################################################
#------------------------- RNAseq + Methyl -------------------------#
gm <- g %du% m
testStatistic <- c("DESeq2", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)")
x=list(rnaseq, imputed_methyl)

fit.iDRWPClass(x=x, y=y, globalGraph=gm, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_17_GM", prob = 0, Gamma = 0, pranking = "t-test", mode = "GM", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GM_18_17 <- fit.classification(y=y, samples = samples, id = "result18_17_GM", datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 50)

save(res_pa_GM_18_17, file=file.path('data/model/res_pa_GM_18_17.RData'))

summary(res_pa_GM_18_17)

write.SigFeatures(res_fit=res_pa_GM_18_17, id = "result18_17_GM", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
gmr <- g %du% m %du% r
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_17_GMR", prob = 0, Gamma = 0, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_18_17 <- fit.classification(y=y, samples = samples, id = "result18_17_GMR", datapath = datapath, respath = respath, 
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_18_17, file=file.path('data/model/res_pa_GMR_18_17.RData'))

summary(res_pa_GMR_18_17)

write.SigFeatures(res_fit=res_pa_GMR_18_17, id = "result18_17_GMR", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_17_GMR_d", prob = 0, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_17 <- fit.classification(y=y, samples = samples, id = "result18_17_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_17, file=file.path('data/model/res_pa_GMR_d_18_17.RData'))

summary(res_pa_GMR_d_18_17)

write.SigFeatures(res_fit=res_pa_GMR_d_18_17, id = "result18_17_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_17_GMP", prob = 0, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_17 <- fit.classification(y=y, samples = samples, id = "result18_17_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_17, file=file.path('data/model/res_pa_GMP_18_17.RData'))

summary(res_pa_GMP_18_17)

write.SigFeatures(res_fit=res_pa_GMP_18_17, id = "result18_17_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GM_18_17, res_pa_GMR_18_17, res_pa_GMR_d_18_17, res_pa_GMP_18_17))


################################################### Result18_18: prob = 0, Gamma = 0.2  #################################################

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_18_GMR_d", prob = 0, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_18 <- fit.classification(y=y, samples = samples, id = "result18_18_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_18, file=file.path('data/model/res_pa_GMR_d_18_18.RData'))

summary(res_pa_GMR_d_18_18)

write.SigFeatures(res_fit=res_pa_GMR_d_18_18, id = "result18_18_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_18_GMP", prob = 0, Gamma = 0.2, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_18 <- fit.classification(y=y, samples = samples, id = "result18_18_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_18, file=file.path('data/model/res_pa_GMP_18_18.RData'))

summary(res_pa_GMP_18_18)

write.SigFeatures(res_fit=res_pa_GMP_18_18, id = "result18_18_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_18_18, res_pa_GMP_18_18))


################################################### Result18_19: prob = 0, Gamma = 0.4  #################################################

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_19_GMR_d", prob = 0, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_19 <- fit.classification(y=y, samples = samples, id = "result18_19_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_19, file=file.path('data/model/res_pa_GMR_d_18_19.RData'))

summary(res_pa_GMR_d_18_19)

write.SigFeatures(res_fit=res_pa_GMR_d_18_19, id = "result18_19_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_19_GMP", prob = 0, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_19 <- fit.classification(y=y, samples = samples, id = "result18_19_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_19, file=file.path('data/model/res_pa_GMP_18_19.RData'))

summary(res_pa_GMP_18_19)

write.SigFeatures(res_fit=res_pa_GMP_18_19, id = "result18_19_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_18_19, res_pa_GMP_18_19))


################################################### Result18_20: prob = 0, Gamma = 0.6  #################################################

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_20_GMR_d", prob = 0, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_20 <- fit.classification(y=y, samples = samples, id = "result18_20_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_20, file=file.path('data/model/res_pa_GMR_d_18_20.RData'))

summary(res_pa_GMR_d_18_20)

write.SigFeatures(res_fit=res_pa_GMR_d_18_20, id = "result18_20_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_20_GMP", prob = 0, Gamma = 0.6, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_20 <- fit.classification(y=y, samples = samples, id = "result18_20_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_20, file=file.path('data/model/res_pa_GMP_18_20.RData'))

summary(res_pa_GMP_18_20)

write.SigFeatures(res_fit=res_pa_GMP_18_20, id = "result18_20_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_18_20, res_pa_GMP_18_20))


################################################### Result18_21: prob = 0, Gamma = 0.8  #################################################

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_21_GMR_d", prob = 0, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_21 <- fit.classification(y=y, samples = samples, id = "result18_21_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_21, file=file.path('data/model/res_pa_GMR_d_18_21.RData'))

summary(res_pa_GMR_d_18_21)

write.SigFeatures(res_fit=res_pa_GMR_d_18_21, id = "result18_21_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_21_GMP", prob = 0, Gamma = 0.8, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_21 <- fit.classification(y=y, samples = samples, id = "result18_21_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_21, file=file.path('data/model/res_pa_GMP_18_21.RData'))

summary(res_pa_GMP_18_21)

write.SigFeatures(res_fit=res_pa_GMP_18_21, id = "result18_21_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_18_21, res_pa_GMP_18_21))


################################################### Result18_22: prob = 0.2, Gamma = 0  #################################################

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_22_GMR_d", prob = 0.2, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_22 <- fit.classification(y=y, samples = samples, id = "result18_22_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_22, file=file.path('data/model/res_pa_GMR_d_18_22.RData'))

summary(res_pa_GMR_d_18_22)

write.SigFeatures(res_fit=res_pa_GMR_d_18_22, id = "result18_22_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_22_GMP", prob = 0.2, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_22 <- fit.classification(y=y, samples = samples, id = "result18_22_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_22, file=file.path('data/model/res_pa_GMP_18_22.RData'))

summary(res_pa_GMP_18_22)

write.SigFeatures(res_fit=res_pa_GMP_18_22, id = "result18_22_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_18_22, res_pa_GMP_18_22))


################################################### Result18_23: prob = 0.4, Gamma = 0  #################################################

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_23_GMR_d", prob = 0.4, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_23 <- fit.classification(y=y, samples = samples, id = "result18_23_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_23, file=file.path('data/model/res_pa_GMR_d_18_23.RData'))

summary(res_pa_GMR_d_18_23)

write.SigFeatures(res_fit=res_pa_GMR_d_18_23, id = "result18_23_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_23_GMP", prob = 0.4, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_23 <- fit.classification(y=y, samples = samples, id = "result18_23_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_23, file=file.path('data/model/res_pa_GMP_18_23.RData'))

summary(res_pa_GMP_18_23)

write.SigFeatures(res_fit=res_pa_GMP_18_23, id = "result18_23_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_18_23, res_pa_GMP_18_23))



################################################### Result18_24: prob = 0.6, Gamma = 0  #################################################

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_24_GMR_d", prob = 0.6, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_24 <- fit.classification(y=y, samples = samples, id = "result18_24_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_24, file=file.path('data/model/res_pa_GMR_d_18_24.RData'))

summary(res_pa_GMR_d_18_24)

write.SigFeatures(res_fit=res_pa_GMR_d_18_24, id = "result18_24_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_24_GMP", prob = 0.6, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_24 <- fit.classification(y=y, samples = samples, id = "result18_24_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_24, file=file.path('data/model/res_pa_GMP_18_24.RData'))

summary(res_pa_GMP_18_24)

write.SigFeatures(res_fit=res_pa_GMP_18_24, id = "result18_24_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_18_24, res_pa_GMP_18_24))


################################################### Result18_25: prob = 0.8, Gamma = 0  #################################################

#------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
gmr <- list(g, m, r)
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_25_GMR_d", prob = 0.8, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMR_d_18_25 <- fit.classification(y=y, samples = samples, id = "result18_25_GMR_d", datapath = datapath, respath = respath,
                                         profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                         nFolds = 5, numTops=50, iter = 50)

save(res_pa_GMR_d_18_25, file=file.path('data/model/res_pa_GMR_d_18_25.RData'))

summary(res_pa_GMR_d_18_25)

write.SigFeatures(res_fit=res_pa_GMR_d_18_25, id = "result18_25_GMR_d", profile_name=profile_name, method="DRW", respath=respath)


#------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa)

fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
               datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
               id = "result18_25_GMP", prob = 0.8, Gamma = 0, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)

res_pa_GMP_18_25 <- fit.classification(y=y, samples = samples, id = "result18_25_GMP", datapath = datapath, respath = respath,
                                       profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                                       nFolds = 5, numTops=50, iter = 50)


save(res_pa_GMP_18_25, file=file.path('data/model/res_pa_GMP_18_25.RData'))

summary(res_pa_GMP_18_25)

write.SigFeatures(res_fit=res_pa_GMP_18_25, id = "result18_25_GMP", profile_name=profile_name, method="DRW", respath=respath)

res_models <- c(res_models, list(res_pa_GMR_d_18_25, res_pa_GMP_18_25))


#########################################################################################################################################
# plot



for(i in 1:length(res_models)){
  if(i %% 4 == 1){
    res_gm <- c(res_gm, list(res_models[[i]]))
  }
  if(i %% 4 == 2){
    res_gmr <- c(res_gmr, list(res_models[[i]]))
  }
  if(i %% 4 == 3){
    res_gmr_d <- c(res_gmr_d, list(res_models[[i]]))
  }
  if(i %% 4 == 0){
    res_gmp <- c(res_gmp, list(res_models[[i]]))
  }
}

xlabs <- c("[g=0.2]", "[g=0.4]", "[g=0.6]", "[g=0.8]")

# Plot for GM models
title <- c("Result 18_GM")
#res_gm <- res_gm[1:4]
perf_min <- min(sapply(X = res_gm, FUN = function(x){mean(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gm, FUN = function(x){mean(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gm, perf_min = perf_min-0.02, perf_max = perf_max+0.02)

# Plot for GMR models
title <- c("Result 18_GMR")
res_gmr <- res_gmr[1:4]
perf_min <- min(sapply(X = res_gmr, FUN = function(x){mean(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr, FUN = function(x){mean(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr, perf_min = perf_min-0.02, perf_max = perf_max+0.02)

xlabs <- c("[p=0.2,g=0.2]", "[p=0.2,g=0.4]", "[p=0.2,g=0.6]", "[p=0.2,g=0.8]", "[p=0.4,g=0.2]", "[p=0.4,g=0.4]", "[p=0.4,g=0.6]", "[p=0.4,g=0.8]", "[p=0.6,g=0.2]", "[p=0.6,g=0.4]", "[p=0.6,g=0.6]", "[p=0.6,g=0.8]", "[p=0.8,g=0.2]", "[p=0.8,g=0.4]", "[p=0.8,g=0.6]", "[p=0.8,g=0.8]")

# Plot for GMR_d models
title <- c("Result 18_GMR_d")
perf_min <- min(sapply(X = res_gmr_d, FUN = function(x){mean(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_d, FUN = function(x){mean(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr_d, perf_min = perf_min-0.02, perf_max = perf_max+0.02)

# Plot for GMP models
title <- c("Result 18_GMP")
perf_min <- min(sapply(X = res_gmp, FUN = function(x){mean(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmp, FUN = function(x){mean(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmp, perf_min = perf_min-0.02, perf_max = perf_max+0.02)
