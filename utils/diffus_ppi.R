diffus_ppi <- function(datapath, gene_weight = gene_weight, ppiGraph = p){
  
  rppapath <- file.path(gdacpath, 'mean_imputed_rppa.csv')
  ppi <- p
  
  # assign rppa mean expression to ppi node(intersect with rppa and ppi gene)
  # The number of common gene : 183(RPPA & PPI)
  # RPPA : 188
  # PPI : 24129
  rppa <- read.csv(rppapath, header=T, row.names = 1)
  common_gene <- intersect(rownames(rppa), substring(V(ppi)$name,2))
  rppa <- rppa[common_gene,] 
  
  row.names(rppa) <- paste("p",rownames(rppa),sep = "")
  
  ########################################################################################
  # read clinical information of BRCA patients
  # RPPA & clinical : 937
  # year :3 => RPPA & clinical : 825
  clinical <- read.csv(file.path(gdacpath, 'brca_clinical.csv'), header=T, row.names=1)
  year <- 3
  
  group_cut <- 365*year
  
  survival_all <- clinical$survival
  samples_all <- rownames(clinical)
  
  samples_r <- colnames(rppa)
  samples_r <- gsub('\\_', '.', samples_r)
  samples_r <- substring(samples_r, 1,15)
  colnames(rppa) <- samples_r
  
  clinical <- clinical[samples_r,]
  
  # remove samples whose survival days were not recorded (NA) or wrongly so as negative values
  remove_samples <- which(is.na(clinical$survival))
  remove_samples <- c(remove_samples, which(clinical$survival<0))
  
  # remove samples whose survival days < 3 years and reported as "living"
  remove_samples <- c(remove_samples, which(clinical$survival<=group_cut & clinical$vital_status==1))
  
  samples_r <- samples_r[-remove_samples]
  clinical <- clinical[-remove_samples,]
  
  rppa <- rppa[,samples_r]
  
  write.csv(rppa, file.path(gdacpath,'rppa_POI.csv'))
  
  # split into 353 good / 472 poor group
  good_samples <- which(clinical$survival>group_cut)
  poor_samples <- which(clinical$survival<=group_cut)
  
  y <- list(good_samples, poor_samples)
  
  # tscore calculation(get W0)
  l_rppa <- list(rppa)
  x_norm <- list(0)
  x_stats <- list(0)
  gene_weight <- list(0)
  DEBUG <- FALSE
  for(i in 1:length(l_rppa)) {
    # normalize gene profile
    x_norm[[i]] <- get.geneprofile.norm(l_rppa[[i]])
    
    # statistics for genes
    x_stats[[i]] <- get.genes.stats(x=l_rppa[[i]], x_norm=x_norm[[i]], y=y,
                                    DEBUG=DEBUG, testStatistic=c('t-test'), pname=c('rppa'), datapath=datapath)
    # initialize gene weights
    geneWeight <- -log(x_stats[[i]][,2]+2.2e-16)
    geneWeight[which(is.na(geneWeight))] <- 0
    gene_weight[[i]] <- (geneWeight - min(geneWeight)) / (max(geneWeight) - min(geneWeight))
    
  }
  
  # assign initial weights to the PPI graph
  W0 <- getW0(gene_weight, ppi)
  
  adj <- as.matrix(get.adjacency(ppi))
  
  # PPI graph diffusion using Markov Random Walk
  dfs_ppi <- random.walk(W0, adj, r=0.5, niter=5000, thresh=1e-05)
  print('PPI diffusion complete..')
  #save(dfs_ppi, file=file.path(datapath, paste(c(name,"rda"), collapse='.')))
  
  return(dfs_ppi[[1]]) # return p.inf the stationary distribution of W0
}
