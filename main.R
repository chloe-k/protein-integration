#source("https://bioconductor.org/biocLite.R")

library(KEGGgraph)
library(igraph)
library(ggplot2)
library(annotate)
library(org.Hs.eg.db)
library(diffusr)
library(Matrix)

sapply(file.path("utils",list.files("utils", pattern="*.R")),source)

# make directory
datapath <- file.path('data')
if(!dir.exists(datapath)) dir.create(datapath)

gdacpath <- file.path(datapath, 'BRCA_GDAC')

respath <- file.path('result')
if(!dir.exists(respath)) dir.create(respath)

# read rda data
graphpath <- file.path(datapath,'directGraph.rda')
if(!file.exists(graphpath)) {
  print('directGraph and pathSet do not exist, now creating directGraph and pathSet is start')
  cons_graph(datapath)
}
load(file.path(graphpath))
load(file.path(datapath, 'pathSet.rda'))

# read RNAseq, Methylation data, RPPA data
data_all_path <- file.path(datapath, "data.RData")
if(!file.exists(data_all_path)) read_data(year, datapath)
load(data_all_path)

ppipath <- file.path(datapath, 'ppiGraph.rda')
dppipath <- file.path(datapath, 'DppiGraph.rda')
if(!file.exists(ppipath)) {
  print('ppiGraph does not exist, now creating ppiGraph start')
  cons_ppi(datapath, gdacpath)
}
load(file.path(ppipath))
load(file.path(dppipath))

# dfsppi <- 'diffus_ppi_8rel'
# dfsppipath <- file.path(datapath, paste(c(dfsppi,'rda'), collapse = '.'))
# if(!file.exists(dfsppipath)){
#   print('diffused PPI with RPPA')
#   diffus_ppi(datapath, gdacpath, p, dfsppi)
# }
# load(file.path(dfsppipath))

# global directed pathway graph provided in DRWPClass
#globGpath <- file.path(datapath, paste(c('dfsppi','rda'), collapse = '.'))

load(file.path(datapath, "directGraph.rda"))

# directed pathway graph provided in DRWPClass
g <- directGraph 
V(g)$name <- paste("g",V(g)$name,sep="")

m <- directGraph
V(m)$name <-paste("m",V(m)$name,sep="")

p <- DppiGraph
V(p)$name <-paste("p",V(p)$name,sep="")

y=list(good_samples, poor_samples)
#---------------DRW---------------#
# RNAseq pathway profile


# Methylation pathway profile


# RPPA pathway profile


#---------------iDRW---------------#
# concat directed pathway graphs within each profile
# RNAseq + Methyl
# gm <- g %du% m
# testStatistic <- c("DESeq2", "t-test")
# profile_name <- c("rna", "meth")
# x=list(rnaseq, imputed_methyl) 
# 
# res_rna_meth <- fit.iDRWPClass(x=x, y=y, globalGraph=gm,
#                                testStatistic= testStatistic, profile_name = profile_name,
#                                datapath = datapath, pathSet=pathSet,
#                                method = "DRW", samples = samples, pranking = "t-test",
#                                iter = 10, AntiCorr=FALSE, DEBUG=TRUE)
# 

# RNAseq + Mehtyl + RPPA
#gmp <- g %du% m %du% p
testStatistic <- c("DESeq2", "t-test", "t-test")
profile_name <- c("rna", "meth", "rppa")
gmp <- list(g, m, p)
x=list(rnaseq, imputed_methyl, rppa) 

# iDRW : RNA-seq + methylation + RPPA profiles (all overlapping genes)
res_pa_GMP <- fit.iDRWPClass(x=x, y=y, globalGraph=gmp,
                                    testStatistic= testStatistic, profile_name = profile_name,
                                    datapath = datapath, respath = respath, pathSet=pathSet,
                                    method = "DRW", samples = samples, pranking = "t-test", mode = "GMP",
                                    iter = 10, AntiCorr=FALSE, DEBUG=FALSE)

save(res_pa_GMP, file=file.path('data/model/res_pa_GMP.RData'))

summary(res_pa_GMP)
print(res_pa_GMP$results)
print(res_pa_GMP$resample$Accuracy)

write.SigFeatures(res_fit=res_pa_GMP, profile_name=profile_name, method="DRW", respath=respath)
