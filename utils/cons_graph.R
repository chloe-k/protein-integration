library(KEGGgraph)
library(igraph)
library(ggplot2)


##################################################################################
##################################################################################
#############Constructing Directed Gene-gene GRAPH################################

# This code implements constructing global directed gene-gene graph from all kegg pathway(525 pathways)
#525 pathways -> 327 pathways(It is the number of existing kgml)
kegg_list <- scan("../data/kegg_pathway", what="", sep = "\n")
kgml_list <- list()

for(pathid in kegg_list){
  tryCatch(
    {
      tmp <- paste('../data/KEGG_DB/hsa',pathid,sep="")
      tmp <- paste(tmp,".xml", sep="")
      retrieveKGML(pathwayid = pathid, organism = "hsa", destfile = tmp, method = "wget")
      kgml_list <- append(kgml_list, tmp)
      cat(pathid,'is downloaded!...................')
      
    },
    error=function(e){
      cat(pathid,'is not exist')
      message(e)
    }
    
  )
}


#kgml_list <- paste('../data/KEGG_DB/',exist_list,sep="")
graphs <- list()

for(kgml in kgml_list){
  pathwayG <- parseKGML2Graph(kgml)
  graphs <- append(graphs, pathwayG)
}



merged_G <- mergeGraphs(graphs)
G <- igraph.from.graphNEL(merged_G)


##################################################################################
##################################################################################
#########################Constructing PPI GRAPH###################################
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
ppi_df <- as.data.frame(ppi.degrees)
hist(ppi_df$ppi.degrees, breaks = 1000, xlab = "Degree", main = "Degree-Frequency histogram in PPI")
ppi_log10 <- log10(ppi_df$ppi.degrees)
hist(ppi_log10, breaks = 1000, xlab = "Degree", main = "Degree-Frequency log10-histogram in PPI")
#hub_gene <- ppi.degrees[ppi.degrees>4000]
#barplot(hub_gene, names.arg = ppi.degrees, cex.names=0.7)
chebi <- grep("CHEBI", rownames(ppi_df))
ppi_without_chebi <- ppi_df[-chebi,1]
hist(ppi_without_chebi, breaks = 1000, xlab = "Degree", main = "Degree-Frequency histogram in PPI without CHEBI")
ppi_without_chebi_log10 <- log10(ppi_without_chebi)
hist(ppi_without_chebi_log10, breaks = 1000, xlab = "Degree", main = "Degree-Frequency log10-histogram in PPI without CHEBI")

#intersection of ppi and rppa : 186
ppi_gene <- V(ppi)
rppa_gene <- rownames(rppa)
common_gene <- intersect(names(ppi_gene),rppa_gene)

