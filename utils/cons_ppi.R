cons_ppi <- function(datapath, gdacpath, rppa){
  
  # PathwayCommons9 is updated in 2017/06/29
  # read data
  ppi_profile <- read.csv(file.path(datapath, 'PathwayCommons9.All.hgnc.sif'), header=F, sep="\t", stringsAsFactors = F)
  #rppa <- read.csv(file.path(gdacpath, 'mean_imputed_rppa.csv'), header=T, row.names=1, sep=",", stringsAsFactors = F)
  gene_name_id_map <- read.csv(file.path(gdacpath, 'gene_name_id_map'), skip=29,header=F, row.names=1, stringsAsFactors = F)
  
  # remove chemical compound(CHEBI) interaction in ppi
  chem <- list()
  chem <- which(ppi_profile[2] == 'consumption-controlled-by')
  chem <- c(chem, which(ppi_profile[2] == 'controls-production-of'))
  chem <- c(chem, which(ppi_profile[2] == 'controls-transport-of-chemical'))
  chem <- c(chem, which(ppi_profile[2] == 'chemical-affects'))
  chem <- c(chem, which(ppi_profile[2] == 'reacts-with'))
  chem <- c(chem, grep('CHEBI', ppi_profile[[1]]), value=FALSE)
  chem <- c(chem, grep('CHEBI', ppi_profile[[3]]), value=FALSE)
  ppi <- ppi_profile[-chem,]
  
  # convert gene symbol to gene id(Entrez ID)(applied after result 9)
  ppi_A <- gene_name_id_map[ppi$V1,]
  unmapped_A <- ppi$V1[which(is.na(ppi_A))]
  ppi_B <- gene_name_id_map[ppi$V3,]
  unmapped_B <- ppi$V3[which(is.na(ppi_B))]
  unmapped_gene_name <- union(unmapped_A, unmapped_B) # unmapped gene name : 3154
  
  ppi$V1 <- gene_name_id_map[ppi$V1,]
  ppi$V3 <- gene_name_id_map[ppi$V3,]
  cat('The number of unmapped gene is ')
  cat(length(unmapped_gene_name))
  cat('\n')
  
  # remove edges(row) which contain unmapped gene name to gene id
  ppi <- na.omit(ppi) # The number of edges 987778  -> 890697 
  
  # bidirected interaction in ppi
  bidir <- list()
  bidir <- which(ppi[2] == 'in-complex-with')
  bidir <- c(bidir, which(ppi[2] == 'interacts-with'))
  bidir <- c(bidir, which(ppi[2] == 'neighbor-of'))
  
  # make reverse direction
  bidir_prot <- ppi[bidir,]
  revdir <- data.frame(matrix(vector(), length(bidir_prot[[1]]), 3, dimnames=list(c(), c("V1", "V2", "V3"))))
  revdir$V1 <- bidir_prot$V3
  revdir$V3 <- bidir_prot$V1
  dppi <- rbind(ppi, revdir)
  
  # make ppi igraph
  ppi[2] <- NULL
  dppi[2] <- NULL
  genepair <- ppi[!duplicated(ppi), ]
  dgenepair <- dppi[!duplicated(dppi),]
  
  #write.table(x=dgenepair, file=file.path(datapath,'dPPI_PathwayCommons9.sif'), sep="\t")
  write.table(x=dgenepair, file=file.path(datapath,'dPPI_PathwayCommons9(Entrez).sif'), sep="\t")
  
  # undirected ppi
  ppiGraph <- graph.data.frame(genepair, directed=FALSE)
  
  # directed ppi
  DppiGraph <- graph.data.frame(dgenepair, directed=TRUE)
  
  # extract subgraph whose nodes are corresponding to RPPA protein(node) and its neighbors.
  ppi_rppa_gene <- intersect(V(ppiGraph)$name, substring(rownames(rppa), 2))
  ppiGraph <- induced_subgraph(ppiGraph, match(ppi_rppa_gene, V(ppiGraph)$name))
  
  Dppi_rppa_gene <- intersect(V(DppiGraph)$name, substring(rownames(rppa), 2))
  DppiGraph <- induced_subgraph(DppiGraph, match(Dppi_rppa_gene, V(DppiGraph)$name))

  #---------------result 1 ~ 8----------------------#  
  # ppi genes : 20057 
  # ppi edges : 904706  (regardless of interaction type) 
  # gsize(ppiGraph)
  
  # ppi genes : 20057 
  # ppi edges : 1360850  (regardless of interaction type) 
  # gsize(DppiGraph)
  #---------------result 1 ~ 8----------------------#
  
  #---------------result 9 ~ ----------------------#
  # ppi genes : 
  # ppi edges : 
  gsize(ppiGraph)
  
  # ppi genes : 
  # ppi edges :  
  gsize(DppiGraph)
  #---------------result 9 ~ ----------------------#
  
  # save(ppiGraph, file=file.path(datapath, paste(c("ppiGraph","rda"), collapse='.')))
  # save(DppiGraph, file=file.path(datapath, paste(c("DppiGraph","rda"), collapse='.')))
  save(ppiGraph, file=file.path(datapath, paste(c("ppiGraph(Entrez)","rda"), collapse='.')))
  save(DppiGraph, file=file.path(datapath, paste(c("DppiGraph(Entrez)","rda"), collapse='.')))
  
  ########################################################################################  
  #plot the degree of the igraph(X:genes, Y:edge)
  #number of the Node degree > 2000 => 13
  ########LIST##########
  #HNF4A  APP   JUN   MAX   MYC   SP1  TP53   SREBF1  TCF3  LEF1  MAZ  FOXO4  NOG 
  ######################
  # ppiGraph.degrees <- degree(ppiGraph)
  # ppi_df <- as.data.frame(ppiGraph.degrees)
  # hist(ppi_df$ppiGraph.degrees, breaks = 1000, xlab = "Degree", main = "Degree-Frequency histogram in PPI")
  # ppi_log10 <- log10(ppi_df$ppiGraph.degrees)
  # hist(ppi_log10, breaks = 1000, xlab = "Degree", main = "Degree-Frequency log10-histogram in PPI")
  # hub_gene <- ppiGraph.degrees[ppiGraph.degrees>2000]
  # #barplot(hub_gene, names.arg = ppiGraph.degrees, cex.names=0.7)
  # 
  # #intersection of ppi and rppa : 186
  # ppi_gene <- V(ppiGraph)
  # rppa_gene <- rownames(rppa)
  # common_gene <- intersect(names(ppi_gene),rppa_gene)
  
}

