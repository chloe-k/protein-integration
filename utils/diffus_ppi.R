diffus_ppi <- function(datapath, gene_weight = gene_weight, ppi, prob){

  # assign initial weights to the PPI graph
  W0 <- getW0(gene_weight, ppi)
  
  adj <- as.matrix(get.adjacency(ppi))
  # adj <- as.matrix(get.adjacency(ppi, attr = "combined_score"))/1000

  ##################### PPI graph diffusion using Markov Random Walk #########################

  # dfs_ppi <- random.walk(W0, adj, r=prob, do.analytical=TRUE)
  # names(dfs_ppi[[1]]) <- names(W0)
  # print('PPI diffusion complete..')
  # 
  # return(dfs_ppi[[1]]) # return p.inf the stationary distribution of W0


  # ##################### PPI graph diffusion using Heat Diffusion Process on Laplacian Matrix #########################
  dfs_ppi <- heat.diffusion(h0 = W0, graph = adj, t = 0.5)
  names(dfs_ppi) <- names(W0)
  print('PPI diffusion complete..')
  print(dfs_ppi)
  return(dfs_ppi)
}

# weighted_kegg <- function(r, p){
#   common_edge <- 
# }
