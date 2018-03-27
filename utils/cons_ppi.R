library(igraph)
library(ggplot2)

#pathway : 1546602 * 2
#rppa : 188 * 376
pathway <- read.csv(file="~/protein-integration/data/PathwayCommons9.All.hgnc.sif", sep="\t", header = FALSE)
rppa <- read.csv(file="~/protein-integration/data/mean_imputed_rppa.csv", header=TRUE, row.names = 1)
substring(rownames(rppa),1,1)
row.names(rppa) <- substring(rownames(rppa),2)

pathway[2] <- NULL
genepair <- unique(pathway)
adjmtx <- get.adjacency(graph.edgelist(as.matrix(genepair), directed=FALSE))
ppi <- graph_from_adjacency_matrix(adjmtx, mode = "undirected")

#ppi genes : 31693
#ppi edges : 1461050 (regardless of interaction type) 
gsize(ppi)

#plot the degree of the igraph(X:gene, Y:edge)
#number of the Node degree > 4000 => 25
########LIST##########
#"APP"         "MYC"         "SP1"         "TCF3"        "CHEBI:16469" "CHEBI:23965" "CHEBI:15367"
#"CHEBI:26536" "CHEBI:33364" "NOG"         "CHEBI:39867" "CHEBI:46195" "CHEBI:23414" "CHEBI:2504" 
#"CHEBI:27899" "CHEBI:29865" "CHEBI:30563" "CHEBI:31440" "CHEBI:4031"  "CHEBI:46024" "CHEBI:4667" 
#"CHEBI:60654" "CHEBI:78510" "CHEBI:91108" "CHEBI:9925" 
######################
ppi.degrees <- degree(ppi)
hub_gene <- ppi.degrees[ppi.degrees>4000]
barplot(hub_gene, names.arg = names(hub_gene), cex.names=0.7)

#intersection of ppi and rppa : 186
ppi_gene <- V(ppi)
rppa_gene <- rownames(rppa)
common_gene <- intersect(names(ppi_gene),rppa_gene)
