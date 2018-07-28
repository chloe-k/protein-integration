preproc_STRING_PPI <- function(datapath){
  
  # Convert Ensembl ID to Entrez gene ID(NCBI gene ID)
  STR_ppi <- read.delim(file = 'data/PPI_DB/9606.protein.links.v10.5.txt', sep = ' ', header = T, stringsAsFactors = F)
  ensembl_mapping <- read.delim(file = 'data/PPI_DB/ensemblTogeneId.tsv', header = T, stringsAsFactors = F, row.names = 2)
  # STR_act <- read.delim(file = 'data/PPI_DB/9606.protein.actions.v10.5.txt', header = T, stringsAsFactors = F)
  # colnames(STR_act) <- c("protein1", "protein2", "mode", "action", "is_directional", "a_is_acting", "score")
  # STR_ppi <- merge(STR, STR_act)
  # STR_ppi <- STR_ppi[,1:3]
  
  # Check unmapped ensembl ID
  protein_A <- ensembl_mapping[STR_ppi$protein1,]
  unmapped_A <- STR_ppi$protein1[which(is.na(protein_A))]
  protein_B <- ensembl_mapping[STR_ppi$protein2,]
  unmapped_B <- STR_ppi$protein2[which(is.na(protein_B))]
  unmapped_gene_name <- union(unmapped_A, unmapped_B) # unmapped ensembl id : 2134
  unmapped_gene_ind <- union(which(is.na(protein_A)), which(is.na(protein_B)))
  
  # Remove rows which contain unmapped ensembl ID
  ppi <- STR_ppi[-unmapped_gene_ind,]
  ppi$protein1 <- ensembl_mapping[ppi$protein1,]
  ppi$protein2 <- ensembl_mapping[ppi$protein2,]
  
  # Construct directed ppi(unweighted)
  DppiGraph <- graph.data.frame(ppi[,1:2], directed=TRUE)
  
  save(DppiGraph, file=file.path(datapath, paste(c("DppiGraph_str","rda"), collapse='.')))
  
  # Construct directed ppi(weighted)
  DppiGraph <- graph.data.frame(ppi, directed=TRUE)

  save(DppiGraph, file=file.path(datapath, paste(c("DppiGraph_W_str","rda"), collapse='.')))
  
  # Construct weighted directGraph(25792 edges have weight)
  # df_directGraph <- as.data.frame(x = get.edgelist(directGraph))
  # colnames(df_directGraph) <- c("protein1", "protein2")
  # merged_directGraph <- merge(df_directGraph, ppi, all.x = TRUE)
  # merged_directGraph$combined_score[is.na(merged_directGraph$combined_score)] <- 0
  pathgraph <- get.data.frame(x = directGraph, what = "both")
  
  dppiedge <- get.edgelist(DppiGraph)
  dppiatt <- get.edge.attribute(DppiGraph)
  dppi <- cbind(dppiedge, dppiatt$combined_score)
  dppi_df <- as.data.frame(dppi)
  dppi_rdc <- dppi[-which(!(dppi_df$V1 %in% pathgraph$vertices[[1]] | dppi_df$V2 %in% pathgraph$vertices[[1]])),]
  
  merge_weight <- merge(x = as.data.frame(pathgraph$edges), y = dppi_rdc, all.x = TRUE)

  directGraph_w <- directGraph
  directGraph_w <- set_edge_attr(directGraph_w, "combined_score", )
  
  # directGraph <- graph.data.frame(merged_directGraph, directed = TRUE)
  
  # save(directGraph, file=file.path(datapath, paste(c("directGraph_W_str","rda"), collapse='.')))
  
  ###############################################################################################################
  # For extracting subgraph which contains only specific interaction type
  # These subgraphs have edge weight
  
  # ppi_1 <- ppi[which(ppi$mode == 'binding'),]
  # ppi_2 <- ppi[which(ppi$mode == 'reaction'),]
  # ppi_3 <- ppi[which(ppi$mode == 'catalysis'),]
  # ppi_4 <- ppi[which(ppi$mode == 'inhibition'),]
  # ppi_5 <- ppi[which(ppi$mode == 'activation'),]
  # ppi_6 <- ppi[which(ppi$mode == 'ptmod'),]
  # ppi_7 <- ppi[which(ppi$mode == 'expression'),]
  # 
  # DppiGraph21_1 <- graph.data.frame(ppi_1, directed=TRUE)
  # DppiGraph21_2 <- graph.data.frame(ppi_2, directed=TRUE)
  # DppiGraph21_3 <- graph.data.frame(ppi_3, directed=TRUE)
  # DppiGraph21_4 <- graph.data.frame(ppi_4, directed=TRUE)
  # DppiGraph21_5 <- graph.data.frame(ppi_5, directed=TRUE)
  # DppiGraph21_6 <- graph.data.frame(ppi_6, directed=TRUE)
  # DppiGraph21_7 <- graph.data.frame(ppi_7, directed=TRUE)
  # 
  # a <- decompose(graph = DppiGraph21_1, mode = "weak", min.vertices = 1)
  # b <- decompose(graph = DppiGraph21_2, mode = "weak", min.vertices = 1)
  # c <- decompose(graph = DppiGraph21_3, mode = "weak", min.vertices = 1)
  # d <- decompose(graph = DppiGraph21_4, mode = "weak", min.vertices = 1)
  # e <- decompose(graph = DppiGraph21_5, mode = "weak", min.vertices = 1)
  # f <- decompose(graph = DppiGraph21_6, mode = "weak", min.vertices = 1)
  # g <- decompose(graph = DppiGraph21_7, mode = "weak", min.vertices = 1)
  # 
  # sapply(X = a, FUN = function(x){length(V(x)$name)}) # 75 graph components
  # sapply(X = b, FUN = function(x){length(V(x)$name)}) # 66 graph components
  # sapply(X = c, FUN = function(x){length(V(x)$name)}) # 20 graph components
  # sapply(X = d, FUN = function(x){length(V(x)$name)}) # 28 graph components
  # sapply(X = e, FUN = function(x){length(V(x)$name)}) # 26 graph components
  # sapply(X = f, FUN = function(x){length(V(x)$name)}) # 40 graph components
  # sapply(X = g, FUN = function(x){length(V(x)$name)}) # 74 graph components
  # 
  # # extract the first graph component(the biggest graph component)
  # DppiGraph21_1 <- a[1] # 12961 nodes, 1287082 edges
  # DppiGraph21_2 <- b[1] # 6816 nodes, 1326062 edges
  # DppiGraph21_3 <- c[1] # 6847 nodes, 801674 edges
  # DppiGraph21_4 <- d[1] # 3402 nodes, 122346 edges
  # DppiGraph21_5 <- e[1] # 7680 nodes, 193286 edges
  # DppiGraph21_6 <- f[1] # 3964 nodes, 71862 edges
  # DppiGraph21_7 <- g[1] # 3465 nodes, 22412 edges
  # 
  # sapply(X = DppiGraph21_1, FUN = function(x){length(intersect(V(x)$name, substring(rownames(rppa),2)))}) # 158 rppa proteins are contained
  # sapply(X = DppiGraph21_2, FUN = function(x){length(intersect(V(x)$name, substring(rownames(rppa),2)))}) # 133 rppa proteins are contained
  # sapply(X = DppiGraph21_3, FUN = function(x){length(intersect(V(x)$name, substring(rownames(rppa),2)))}) # 128 rppa proteins are contained
  # sapply(X = DppiGraph21_4, FUN = function(x){length(intersect(V(x)$name, substring(rownames(rppa),2)))}) # 113 rppa proteins are contained
  # sapply(X = DppiGraph21_5, FUN = function(x){length(intersect(V(x)$name, substring(rownames(rppa),2)))}) # 150 rppa proteins are contained
  # sapply(X = DppiGraph21_6, FUN = function(x){length(intersect(V(x)$name, substring(rownames(rppa),2)))}) # 134 rppa proteins are contained
  # sapply(X = DppiGraph21_7, FUN = function(x){length(intersect(V(x)$name, substring(rownames(rppa),2)))}) # 137 rppa proteins are contained
  # 
  # save(DppiGraph21_1, file=file.path(datapath, paste(c("DppiGraph21_1","rda"), collapse='.')))
  # save(DppiGraph21_2, file=file.path(datapath, paste(c("DppiGraph21_2","rda"), collapse='.')))
  # save(DppiGraph21_3, file=file.path(datapath, paste(c("DppiGraph21_3","rda"), collapse='.')))
  # save(DppiGraph21_4, file=file.path(datapath, paste(c("DppiGraph21_4","rda"), collapse='.')))
  # save(DppiGraph21_5, file=file.path(datapath, paste(c("DppiGraph21_5","rda"), collapse='.')))
  # save(DppiGraph21_6, file=file.path(datapath, paste(c("DppiGraph21_6","rda"), collapse='.')))
  # save(DppiGraph21_7, file=file.path(datapath, paste(c("DppiGraph21_7","rda"), collapse='.')))
  ###############################################################################################################
}




