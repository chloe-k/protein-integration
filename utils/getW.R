getW <- function(G, gene_weight, x, datapath, mode, EdgeWeight=FALSE, AntiCorr=FALSE){
  print('getting adjacency matrix starts(getW)..')
  len <- length(gene_weight)
  
  if(len > 1){
    if(!AntiCorr) {
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
      } else if(mode == "GMR_1"){
        # make Global Graph(RNAseq+Mehtyl+RPPA(Pathway Graph))
        intersect_mr <- intersect(substring(names(gene_weight[[2]]),2), substring(names(gene_weight[[3]]),2))
        intersect_gr <- intersect(substring(names(gene_weight[[1]]),2), substring(names(gene_weight[[3]]),2))
        
        # create Global Graph's adjacency matrix
        W <- as.matrix(get.adjacency(G)) 
        
        # add edges (RPPA->Methyl)
        for(i in 1:length(intersect_mr)){
          idx <- which(paste("p", intersect_mr[i], sep="") == rownames(W))
          if(length(idx) > 0){
            W[paste("p",intersect_mr[i],sep=""),paste("m",intersect_mr[i],sep="")] <- 1
          }
        }
        
        # add edges (RPPA->RNAseq)
        for(i in 1:length(intersect_gr)){
          idx <- which(paste("g", intersect_gr[i], sep="") == rownames(W))
          if(length(idx) > 0){
            W[paste("p",intersect_gr[i],sep=""),paste("g",intersect_gr[i],sep="")] <- 1
          }
        }
        
      }else if(mode == "GMR_bidir"){
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
            W[paste("g",intersect_gm[i],sep=""),paste("m",intersect_gm[i],sep="")] <- 1
          }
        }
        
        # add edges (RPPA->RNAseq)
        for(i in 1:length(intersect_gr)){
          idx <- which(paste("g", intersect_gr[i], sep="") == rownames(W))
          if(length(idx) > 0){
            W[paste("p",intersect_gr[i],sep=""),paste("g",intersect_gr[i],sep="")] <- 1
            W[paste("g",intersect_gr[i],sep=""),paste("p",intersect_gr[i],sep="")] <- 1
          }
        }
      }
    } else if(AntiCorr) {
      # assign bi-directional edges to significantly anti-correlated genes between exp & meth
      if(mode == 'GM'){
        intersect_genes <- Reduce(intersect, lapply(gene_weight, function(x) substring(names(x),2))) 
        W <- as.matrix(get.adjacency(G))
        xG <- x[[1]]
        xM <- x[[2]]
        
        for(i in 1:length(intersect_genes)){
          idx=which(paste("g",intersect_genes[i],sep="")==rownames(W))
          if(length(idx)>0) {
            if(cor(t(xG[paste("g",intersect_genes[i],sep=""),]),
                   t(xM[paste("m",intersect_genes[i],sep=""),])) < 0 &
               cor.test(t(xG[paste("g",intersect_genes[i],sep=""),]),
                        t(xM[paste("m",intersect_genes[i],sep=""),]),
                        method = "pearson", alternative = "less")$p.value <= 0.05) {
              # W[paste("g",intersect_genes[i],sep=""),paste("m",intersect_genes[i],sep="")] <- 1
              W[paste("m",intersect_genes[i],sep=""),paste("g",intersect_genes[i],sep="")] <- 1
            }
          }
        }
      } else if(mode == "GMR_2"){
        W <- as.matrix(get.adjacency(G))
        
        xG <- x[[1]]
        xM <- x[[2]]
        xP <- x[[3]]
        
        intersect_gm <- intersect(substring(names(gene_weight[[1]]),2), substring(names(gene_weight[[2]]),2))
        intersect_gp <- intersect(substring(names(gene_weight[[1]]),2), substring(names(gene_weight[[3]]),2)) 
        intersect_mp <- intersect(substring(names(gene_weight[[2]]),2), substring(names(gene_weight[[3]]),2))
        
        # add edge (Methyl -> RNA-seq only anticorrelated)
        for(i in 1:length(intersect_gm)){
          idx=which(paste("m",intersect_gm[i],sep="")==rownames(W))
          if(length(idx)>0) {
            if(cor(t(xG[paste("g",intersect_gm[i],sep=""),]),
                   t(xM[paste("m",intersect_gm[i],sep=""),])) < 0 &
               cor.test(t(xG[paste("g",intersect_gm[i],sep=""),]),
                        t(xM[paste("m",intersect_gm[i],sep=""),]),
                        method = "pearson", alternative = "less")$p.value <= 0.05) {
              W[paste("m",intersect_gm[i],sep=""),paste("g",intersect_gm[i],sep="")] <- 1
            }
          }
        }
        
        # add edges (RPPA->RNAseq)
        for(i in 1:length(intersect_gp)){
          idx <- which(paste("p", intersect_gp[i], sep="") == rownames(W))
          if(length(idx) > 0){
            W[paste("p",intersect_gp[i],sep=""),paste("g",intersect_gp[i],sep="")] <- 1
          }
        }
        
        # add edge (RPPA -> Methyl only anticorrelated)
        for(i in 1:length(intersect_mp)){
          idx=which(paste("p",intersect_mp[i],sep="")==rownames(W))
          if(length(idx)>0) {
            if(cor(t(xP[paste("p",intersect_mp[i],sep=""),]),
                   t(xM[paste("m",intersect_mp[i],sep=""),])) < 0 &
               cor.test(t(xP[paste("p",intersect_mp[i],sep=""),]),
                        t(xM[paste("m",intersect_mp[i],sep=""),]),
                        method = "pearson", alternative = "less")$p.value <= 0.05) {
              W[paste("p",intersect_mp[i],sep=""),paste("m",intersect_mp[i],sep="")] <- 1
            }
          }
        }
      }
      
      else if(mode == "GMR_3" | mode == "GMR_3_d"){
        W <- as.matrix(get.adjacency(G))
        
        xG <- x[[1]]
        xM <- x[[2]]
        xP <- x[[3]]
        
        intersect_gm <- intersect(substring(names(gene_weight[[1]]),2), substring(names(gene_weight[[2]]),2))
        intersect_gp <- intersect(substring(names(gene_weight[[1]]),2), substring(names(gene_weight[[3]]),2)) 
        
        
        # add edge (Methyl -> RNA-seq only anticorrelated)
        for(i in 1:length(intersect_gm)){
          idx=which(paste("m",intersect_gm[i],sep="")==rownames(W))
          if(length(idx)>0) {
            if(cor(t(xG[paste("g",intersect_gm[i],sep=""),]),
                   t(xM[paste("m",intersect_gm[i],sep=""),])) < 0 &
               cor.test(t(xG[paste("g",intersect_gm[i],sep=""),]),
                        t(xM[paste("m",intersect_gm[i],sep=""),]),
                        method = "pearson", alternative = "less")$p.value <= 0.05) {
              W[paste("m",intersect_gm[i],sep=""),paste("g",intersect_gm[i],sep="")] <- 1
            }
          }
        }
        
        # add edges (RPPA->RNAseq)
        for(i in 1:length(intersect_gp)){
          idx <- which(paste("p", intersect_gp[i], sep="") == rownames(W))
          if(length(idx) > 0){
            W[paste("p",intersect_gp[i],sep=""),paste("g",intersect_gp[i],sep="")] <- 1
          }
        }
        
      }
      
    }
    
    
    
  }else{
    W <- as.matrix(get.adjacency(G))
  }
  
  
  print(dim(W)) # number of nodes
  print(sum(W)) # number of edges (adjacency matrix)

  W[is.na(W)] <- 0
  save(W, file=file.path(datapath, paste(c(mode, if(AntiCorr) "anticorr","W","RData"), collapse = '.')))
  print('Adjacency matrix W complete ...')
  
  return(W)
  
}