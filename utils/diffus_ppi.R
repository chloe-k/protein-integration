diffus_ppi <- function(datapath, gene_weight = gene_weight, ppi, prob){

  # assign initial weights to the PPI graph
  W0 <- getW0(gene_weight, ppi)
  
  # adj <- as.matrix(get.adjacency(ppi))
  adj <- as.matrix(get.adjacency(ppi, attr = "combined_score"))/100
  
  ##################### PPI graph diffusion using Markov Random Walk #########################
  # dfs_ppi <- random.walk(W0, adj, r=0.5, niter=5000, thresh=1e-05)
  dfs_ppi <- random.walk(W0, adj, r=prob, niter=5000, thresh=1e-05)
  names(dfs_ppi[[1]]) <- names(W0)
  print('PPI diffusion complete..')
  #save(dfs_ppi, file=file.path(datapath, paste(c(name,"rda"), collapse='.')))

  return(dfs_ppi[[1]]) # return p.inf the stationary distribution of W0


  # ##################### PPI graph diffusion using Heat Diffusion Process on Laplacian Matrix #########################
  # dfs_ppi <- heat.diffusion(W0, adj)
  # names(dfs_ppi) <- names(W0)
  # print('PPI diffusion complete..')
  # return(dfs_ppi)
}
