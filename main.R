source("https://bioconductor.org/biocLite.R")
<<<<<<< HEAD
=======

>>>>>>> 013e1e2f2309762dcabe59dbcde9edb88eca9def
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

# imputation ppi node value by using RPPA profile

rppapath <- file.path(datapath, 'mean_imputed_rppa.csv')

if(!file.exists(rppapath)) {
  print('preprocessed rppa profile does not exist, now preprocessing RPPA profile start')
  preprcs_rppa(gdacpath)
}

ppipath <- file.path(datapath, 'ppiGraph.rda')
dppipath <- file.path(datapath, 'DppiGraph.rda')
if(!file.exists(ppipath)) {
  print('ppiGraph does not exist, now creating ppiGraph start')
  cons_ppi(datapath, gdacpath)
}
load(file.path(ppipath))
load(file.path(dppipath))


g <- directGraph
m <- directGraph
dp <- DppiGraph
p <- ppiGraph

dfsppi <- 'diffus_ppi_8rel'
dfsppipath <- file.path(datapath, paste(c(dfsppi,'rda'), collapse = '.'))
if(!file.exists(dfsppipath)){
  print('diffused PPI with RPPA')
  diffus_ppi(datapath, gdacpath, p, dfsppi)
}
load(file.path(dfsppipath))






