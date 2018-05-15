#CPTAC
#10628*228
cptac <- read.csv(file="~/protein-integration/data/TCGA_Breast_BI_Proteome.itraq.tsv", sep="\t", header = T, row.names = 1)

#PPI
pathway <- read.csv(file="~/protein-integration/data/PathwayCommons9.All.hgnc.sif", sep="\t", header = FALSE)

#preprocessing data
cptac <- cptac[-(1:3),]

pathway[2] <- NULL
genepair <- unique(pathway)
adjmtx <- get.adjacency(graph.edgelist(as.matrix(genepair), directed=FALSE))

#construct ppi
ppi <- graph_from_adjacency_matrix(adjmtx, mode = "undirected")
ppi.degrees <- degree(ppi)

#check intersect of cptac and ppi genes
#10403
length(intersect(names(ppi.degrees), rownames(cptac)))

#Pathway graph & RPPA gene : 197
rppa <- read.csv(file="~/protein-integration/data/rppa_POI.csv", header = T, row.names = 1, stringsAsFactors = F)
rppa_pathSet <- list()

sink("data/pathway_rppa.txt")
cat("Pathway\tnum_RPPA_genes\tRPPA_gene\n")

for(i in 1:length(pathSet)){
  len_rpath <- length(intersect(pathSet[[i]],rownames(rppa)))
  if(len_rpath>0){
    rppa_pathSet[[names(pathSet[i])]] <- unlist(intersect(pathSet[[i]],rownames(rppa)))
    num_rppa <- length(unlist(intersect(pathSet[[i]],rownames(rppa))))
    cat(names(pathSet[i]))
    cat('\t')
    cat(num_rppa)
    cat('\t')
    for(cnt in 1:num_rppa){
      cat(unlist(intersect(pathSet[[i]],rownames(rppa)))[cnt])
      cat('\t')
    }
    cat('\n')
  }
}
sink()

