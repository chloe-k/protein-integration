cons_ppi <- function(datapath, from){
  #library(igraph)
  #library(ggplot2)
  
  #pathway : 1546602 * 2
  #rppa : 188 * 376
  #PathwayCommons9 is updated in 2017/06/29
  
  datapath <- '~/protein-integration/data/'
  #pathway <- read.csv(file="~/protein-integration/data/PathwayCommons9.All.hgnc.sif", sep="\t", header = FALSE)
  #rppa <- read.csv(file="~/protein-integration/data/mean_imputed_rppa.csv", header=TRUE, row.names = 1)
  pathway <- read.csv(file.path(datapath, 'PathwayCommons9.All.hgnc.sif'), header=F, sep="\t")
  rppa <- read.csv(file.path(datapath, 'mean_imputed_rppa.csv'), header=T, row.names=1, sep=",")
  
  #remove the prefix 'r' of row name
  substring(rownames(rppa),1,1)
  row.names(rppa) <- substring(rownames(rppa),2)
  
  #remove chemical compound interaction in ppi
  chem <- which(pathway[2] == 'consumption-controlled-by')
  chem <- c(chem, which(pathway[2] == 'controls-production-of'))
  chem <- c(chem, which(pathway[2] == 'controls-transport-of-chemical'))
  chem <- c(chem, which(pathway[2] == 'chemical-affects'))
  chem <- c(chem, which(pathway[2] == 'reacts-with'))
  pathway <- pathway[-chem,]
  
  pathway[2] <- NULL
  genepair <- unique(pathway)
  adjmtx <- get.adjacency(graph.edgelist(as.matrix(genepair), directed=FALSE))
  ppi <- graph_from_adjacency_matrix(adjmtx, mode = "undirected")
  
  #ppi proteins : 24129 
  #ppi edges : 919192  (regardless of interaction type) 
  gsize(ppi)
  
  #plot the degree of the igraph(X:protein, Y:edge)
  #number of the Node degree > 2000 => 13
  ########LIST##########
  #HNF4A  APP   JUN   MAX   MYC   SP1  TP53   SREBF1  TCF3  LEF1  MAZ  FOXO4  NOG 
  ######################
  ppi.degrees <- degree(ppi)
  ppi_df <- as.data.frame(ppi.degrees)
  hist(ppi_df$ppi.degrees, breaks = 1000, xlab = "Degree", main = "Degree-Frequency histogram in PPI")
  ppi_log10 <- log10(ppi_df$ppi.degrees)
  hist(ppi_log10, breaks = 1000, xlab = "Degree", main = "Degree-Frequency log10-histogram in PPI")
  hub_protein <- ppi.degrees[ppi.degrees>2000]
  barplot(hub_protein, names.arg = ppi.degrees, cex.names=0.7)
  
  #intersection of ppi and rppa : 186
  ppi_gene <- V(ppi)
  rppa_gene <- rownames(rppa)
  common_gene <- intersect(names(ppi_gene),rppa_gene)
  
}

