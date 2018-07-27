getW <- function(datapath, G, gene_weight, mode){
  print('getting adjacency matrix starts(getW)..')
  len <- length(gene_weight)
  if(len > 1){
    if(mode == "GM"){
      # make Global Graph(RNAseq+Mehtyl)
      intersect_gm <- Reduce(intersect, lapply(gene_weight, function(x) substring(names(x),2)))
      
      # create Global Graph's adjacency matrix
      W <- as.matrix(get.adjacency(G)) 
      for(i in 1:length(intersect_gm)){
        idx <- which(paste("g", intersect_gm[i], sep="") == rownames(W))
        if(length(idx) > 0){
          W[paste("m",intersect_gm[i],sep=""),paste("g",intersect_gm[i],sep="")] <- 1
        }
      }
      
      
    } else if(mode == "GR" | mode == "GR_d"){
      # make Global Graph(RNAseq+RPPA)
      intersect_gp <- Reduce(intersect, lapply(gene_weight, function(x) substring(names(x),2)))
      
      # create Global Graph's adjacency matrix
      W <- as.matrix(get.adjacency(G)) 
      for(i in 1:length(intersect_gp)){
        idx <- which(paste("g", intersect_gp[i], sep="") == rownames(W))
        if(length(idx) > 0){
          W[paste("p",intersect_gp[i],sep=""),paste("g",intersect_gp[i],sep="")] <- 1
        }
      }
    } else if(mode == "GMR"| mode == "GMR_d"){
      # make Global Graph(RNAseq+Mehtyl+RPPA(Pathway Graph))
      intersect_gm <- intersect(substring(names(gene_weight[[1]]),2), substring(names(gene_weight[[2]]),2))
      intersect_gr <- intersect(substring(names(gene_weight[[1]]),2), substring(names(gene_weight[[3]]),2))
      
      # create Global Graph's adjacency matrix
      W <- as.matrix(get.adjacency(G)) 
      
      # add edges (Methyl->RNAseq)
      for(i in 1:length(intersect_gm)){
        idx <- which(paste("g", intersect_gm[i], sep="") == rownames(W))
        if(length(idx) > 0){
          W[paste("m",intersect_gm[i],sep=""),paste("g",intersect_gm[i],sep="")] <- 1
        }
      }
      
      # add edges (RPPA->RNAseq)
      for(i in 1:length(intersect_gr)){
        idx <- which(paste("g", intersect_gr[i], sep="") == rownames(W))
        if(length(idx) > 0){
          W[paste("p",intersect_gr[i],sep=""),paste("g",intersect_gr[i],sep="")] <- 1
        }
      }
      
    } else if(mode == "GMP"){
      # make Global Graph(RNAseq+Methyl+RPPA)
      # gene_weight[[1]] : RNAseq(G)
      # gene_weight[[2]] : Methyl(M)
      # gene_weight[[3]] : RPPA(P)
      
      intersect_gm <- intersect(substring(names(gene_weight[[1]]),2), substring(names(gene_weight[[2]]),2))
      intersect_gp <- intersect(substring(names(gene_weight[[1]]),2), substring(names(gene_weight[[3]]),2)) 
      
      # create Global Graph's adjacency matrix
      W <- as.matrix(get.adjacency(G)) 
      
      # add edges (Methyl->RNAseq)
      for(i in 1:length(intersect_gm)){
        idx <- which(paste("g", intersect_gm[i], sep="") == rownames(W))
        if(length(idx) > 0){
          W[paste("m",intersect_gm[i],sep=""),paste("g",intersect_gm[i],sep="")] <- 1 # m->g
        }
      }
      
      # add edges (RPPA->RNAseq)
      for(i in 1:length(intersect_gp)){
        idx <- which(paste("g", intersect_gp[i], sep="") == rownames(W))
        if(length(idx) > 0){
          W[paste("p",intersect_gp[i],sep=""),paste("g",intersect_gp[i],sep="")] <- 1 # p->g
        }
      }
    }
    
  }else{
    W <- as.matrix(get.adjacency(G))
  }
  
  
  print(dim(W)) # number of nodes
  print(sum(W)) # number of edges (adjacency matrix)

  W[is.na(W)] <- 0
  save(W, file=file.path(datapath, paste(c(mode,"W","RData"), collapse = '.')))
  print('Adjacency matrix W complete ...')
  
  return(W)
  
}