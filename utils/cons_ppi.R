cons_ppi <- function(datapath){
  
  # PathwayCommons9 is updated in 2017/06/29
  
  #datapath <- '~/protein-integration/data/'
  pathway_profile <- read.csv(file.path(datapath, 'PathwayCommons9.All.hgnc.sif'), header=F, sep="\t")
  rppa <- read.csv(file.path(datapath, 'mean_imputed_rppa.csv'), header=T, row.names=1, sep=",")
  
  # remove chemical compound(CHEBI) interaction in ppi
  chem <- list()
  chem <- which(pathway_profile[2] == 'consumption-controlled-by')
  chem <- c(chem, which(pathway_profile[2] == 'controls-production-of'))
  chem <- c(chem, which(pathway_profile[2] == 'controls-transport-of-chemical'))
  chem <- c(chem, which(pathway_profile[2] == 'chemical-affects'))
  chem <- c(chem, which(pathway_profile[2] == 'reacts-with'))
  chem <- c(chem, grep('CHEBI', pathway_profile[[1]]), value=FALSE)
  chem <- c(chem, grep('CHEBI', pathway_profile[[3]]), value=FALSE)
  pathway <- pathway_profile[-chem,]
  
  # bidirected interaction in ppi
  bidir <- list()
  bidir <- which(pathway[2] == 'in-complex-with')
  bidir <- c(bidir, which(pathway[2] == 'interacts-with'))
  bidir <- c(bidir, which(pathway[2] == 'neighbor-of'))
  
  # make reverse direction
  bidir_prot <- pathway[bidir,]
  revdir <- data.frame(matrix(vector(), length(bidir_prot[[1]]), 3, dimnames=list(c(), c("V1", "V2", "V3"))))
  revdir$V1 <- bidir_prot$V3
  revdir$V3 <- bidir_prot$V1
  dpathway <- rbind(pathway, revdir)
  
  # make ppi igraph
  pathway[2] <- NULL
  dpathway[2] <- NULL
  genepair <- pathway[!duplicated(pathway), ]
  dgenepair <- dpathway[!duplicated(dpathway),]
  
  # undirected ppi
  ppiGraph <- graph.data.frame(genepair, directed=FALSE)
  
  # directed ppi
  DppiGraph <- graph.data.frame(dgenepair, directed=TRUE)
  
  # ppi genes : 20057 
  # ppi edges : 904706  (regardless of interaction type) 
  gsize(ppiGraph)
  
  # ppi genes : 20057 
  # ppi edges : 1360850  (regardless of interaction type) 
  gsize(DppiGraph)
  
  save(ppiGraph, file=file.path(datapath, paste(c("ppiGraph","rda"), collapse='.')))
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

