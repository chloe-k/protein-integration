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
