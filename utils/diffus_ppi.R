diffus_ppi <- function(datapath, gene_weight = gene_weight, ppi){

  # assign initial weights to the PPI graph
  W0 <- getW0(gene_weight, ppi)
  
  adj <- as.matrix(get.adjacency(ppi))
  adj <- t(adj)
  # PPI graph diffusion using Markov Random Walk
  dfs_ppi <- random.walk(W0, adj, r=0.5, niter=5000, thresh=1e-05)
  names(dfs_ppi[[1]]) <- names(W0)
  print('PPI diffusion complete..')
  #save(dfs_ppi, file=file.path(datapath, paste(c(name,"rda"), collapse='.')))
  
  return(dfs_ppi[[1]]) # return p.inf the stationary distribution of W0
}
