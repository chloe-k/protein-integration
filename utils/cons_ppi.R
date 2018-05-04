cons_ppi <- function(datapath){
  
  #PathwayCommons9 is updated in 2017/06/29
  
  datapath <- '~/protein-integration/data/'
  pathway <- read.csv(file.path(datapath, 'PathwayCommons9.All.hgnc.sif'), header=F, sep="\t")
  rppa <- read.csv(file.path(datapath, 'mean_imputed_rppa.csv'), header=T, row.names=1, sep=",")
  
  #remove the prefix 'r' of row name
  substring(rownames(rppa),1,1)
  row.names(rppa) <- substring(rownames(rppa),2)
  
  #remove chemical compound(CHEBI) interaction in ppi
  chem <- which(pathway[2] == 'consumption-controlled-by')
  chem <- c(chem, which(pathway[2] == 'controls-production-of'))
  chem <- c(chem, which(pathway[2] == 'controls-transport-of-chemical'))
  chem <- c(chem, which(pathway[2] == 'chemical-affects'))
  chem <- c(chem, which(pathway[2] == 'reacts-with'))
  pathway <- pathway[-chem,]
  
  #make ppi igraph
  pathway[2] <- NULL
  genepair <- unique(pathway)
  adjmtx <- get.adjacency(graph.edgelist(as.matrix(genepair), directed=FALSE))
  #ppiGraph <- graph_from_adjacency_matrix(adjmtx, mode = "undirected")
  DppiGraph <- graph_from_adjacency_matrix(adjmtx, mode = "directed")
  
  #ppi genes : 24129 
  #ppi edges : 919192  (regardless of interaction type) 
  #gsize(ppiGraph)
  
  #directed ppi genes : 31693 
  #directed ppi edges : 2922100 
  gsize(DppiGraph)
  
  #save(ppiGraph, file=file.path(datapath, paste(c("ppiGraph","rda"), collapse='.')))
  save(DppiGraph, file=file.path(datapath, paste(c("DppiGraph","rda"), collapse='.')))

########################################################################################  
  #plot the degree of the igraph(X:genes, Y:edge)
  #number of the Node degree > 2000 => 13
  ########LIST##########
  #HNF4A  APP   JUN   MAX   MYC   SP1  TP53   SREBF1  TCF3  LEF1  MAZ  FOXO4  NOG 
  ######################
  ppiGraph.degrees <- degree(ppiGraph)
  ppi_df <- as.data.frame(ppiGraph.degrees)
  hist(ppi_df$ppiGraph.degrees, breaks = 1000, xlab = "Degree", main = "Degree-Frequency histogram in PPI")
  ppi_log10 <- log10(ppi_df$ppiGraph.degrees)
  hist(ppi_log10, breaks = 1000, xlab = "Degree", main = "Degree-Frequency log10-histogram in PPI")
  hub_gene <- ppiGraph.degrees[ppiGraph.degrees>2000]
  #barplot(hub_gene, names.arg = ppiGraph.degrees, cex.names=0.7)
  
  #intersection of ppi and rppa : 186
  ppi_gene <- V(ppiGraph)
  rppa_gene <- rownames(rppa)
  common_gene <- intersect(names(ppi_gene),rppa_gene)
  
}

