library(KEGGgraph)
library(igraph)
library(ggplot2)
library(annotate)
source("https://bioconductor.org/biocLite.R")
library(org.Hs.eg.db)

datapath <- file.path('data')
if(!dir.exists(datapath)) dir.create(datapath)

