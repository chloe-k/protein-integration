fit.iDRWPClass <-
  function(x, y, testStatistic, profile_name, globalGraph = NULL, datapath, respath, pathSet,
           method = "DRW", samples, pranking = "t-test", mode, lim = NULL,
           classifier, nFolds, numTops, id, prob, type_used = NULL,
           iter, Gamma, AntiCorr = FALSE, DEBUG=TRUE) {
    
    if(mode == 'GMR' | mode == 'GM'){
      subId <- paste(c(mode,'_g',Gamma), collapse = '')
    }else if(mode == 'GMR_d'){
      subId <- paste(c(mode,'_p',prob,'_g',Gamma), collapse = '')
    }
    pathAct_param_path <- file.path(datapath, paste(c("pathAct_param_", subId, ".RData"), collapse = ''))
    
    if(!file.exists(pathAct_param_path)){
      cat('Parameters for calculating pathway activity do not exist\n')
      cat(c(pathAct_param_path, '\n'))
      
      x_norm <- list(0)
      x_stats <- list(0)
      gene_weight <- list(0)
      
      for(i in 1:length(x)) {
        # normalize gene profile
        x_norm[[i]] <- get.geneprofile.norm(x[[i]])
        
        # statistics for genes
        x_stats[[i]] <- get.genes.stats(x=x[[i]], x_norm=x_norm[[i]], y=y, 
                                        DEBUG=DEBUG, testStatistic=testStatistic[i], pname=profile_name[i], datapath=datapath)
        # initialize gene weights
        geneWeight <- -log(x_stats[[i]][,2]+2.2e-16)
        geneWeight[which(is.na(geneWeight))] <- 0
        gene_weight[[i]] <- (geneWeight - min(geneWeight)) / (max(geneWeight) - min(geneWeight))
        
      }
      
      if(method == "DRW") {
        # assign initial weights to the pathway graph
        if(mode == "GMP" | mode == "GMR_d"){
          
          # get W0 of G & M 
          gm <- globalGraph[[1]] %du% globalGraph[[2]]
          W0 <- getW0(list(gene_weight[[1]], gene_weight[[2]]), gm)
          
          # get W0 of P
          p_W0 <- diffus_ppi(datapath = datapath, gene_weight = list(gene_weight[[3]]), ppi = globalGraph[[3]], prob = prob)
          
          # concatenate W0 and p_W0
          W0 <- c(W0, p_W0, use.names = TRUE)
          if(DEBUG) cat('Getting W0 is done...')
          
          gmp <- gm %du% globalGraph[[3]]
          
          # get adjacency matrix of the (integrated) gene-gene graph
          wpath <- file.path(datapath, paste(c(mode,"W","RData"), collapse = '.'))
          if(!file.exists(wpath)){
            W = getW(datapath = datapath, G = gmp, gene_weight = gene_weight, mode = mode)
          }
          W = get(load(wpath))
          
        } 
        else{
          W0 <- getW0(gene_weight, globalGraph)
          # W0 <- diffus_ppi(datapath = datapath, gene_weight = gene_weight, ppi = globalGraph, prob = prob)
          if(DEBUG) cat('Getting W0 is done...')
          
          # get adjacency matrix of the (integrated) gene-gene graph
          wpath <- file.path(datapath, paste(c(mode,"W","RData"), collapse = '.'))
          if(!file.exists(wpath)){
            W = getW(datapath = datapath, G = globalGraph, gene_weight = gene_weight, mode = mode)
          }
          W = get(load(wpath))
          
        }
        
        # perform DRW on gene-gene graph
        if(DEBUG) cat('Performing directed random walk...')
        vertexWeight <- DRW(W = W, p0 = W0, gamma = Gamma)
        names(vertexWeight) <- names(W0)
        if(DEBUG) cat('Done\n')	
      } else {
        vertexWeight <- NULL
      }
      
      # reduce list of profile matrices
      x <- Reduce(rbind, x_norm)
      x_stats <- Reduce(rbind, x_stats)
      desc <- c(profile_name, method, if(AntiCorr) "anticorr", "txt")
      
      # save W_inf, x_norm, x_stats for calculating pathway activity score
      save(x, x_stats, vertexWeight, desc, file=pathAct_param_path)
      
    }else{
      cat('Parameters for calculating pathway activity already exist\n')
      cat(c(pathAct_param_path, '\n'))
      load(file = pathAct_param_path)
    }
    
    if(method == "gf") {

      fname_rank = file.path(respath, paste(c("gene_rank", id, testStatistic, desc), collapse = '.'))

      # rank genes by their scores
      stats_feats <- x_stats[ ,1]
      Idx <- sort(stats_feats,decreasing = TRUE,index.return=TRUE)$ix
      stats_feats <- stats_feats[Idx]

      X <- t(x)

      save(stats_feats, X, file=file.path(datapath, paste(c("gene_rank", id, profile_name, method, "RData"), collapse = '.')))

    } else {

      # pathway activity inference method
      # method = DRW / mean / median
      fname_profile = file.path(respath, paste(c("pathway_profile", id, desc), collapse = '.'))
      pApath <- file.path(respath, paste(c("pA", id, profile_name, method, if(AntiCorr) "anticorr", "RData"), collapse = '.'))

        pA <- getPathActivity(x = x, pathSet = pathSet, w = vertexWeight, vertexZP = x_stats, lim=lim,
                              method = method, fname = fname_profile, rows = samples, type_used = type_used)

        save(pA, file=pApath)


      # rank pathway activities
      # ranking = t-test / DA
      fname_rank = file.path(respath, paste(c("pathway_rank", id, pranking, desc), collapse = '.'))

      stats_feats <- rankPathActivity(pathActivity = pA$pathActivity, y = y,
                                      ranking = pranking, fname=fname_rank)

      X <- t(pA$pathActivity)

      save(stats_feats, X, file=file.path(datapath, paste(c("pathway_rank", id, profile_name, method, "RData"), collapse = '.')))
    }

    # write pathway ranking
    write.table(x=matrix(stats_feats, nrow=length(stats_feats), dimnames=list(names(stats_feats),"rank")),
                file=fname_rank, sep="\t", row.names=T, col.names=T)


    cat('fit.iDRWPClass is done\n')
    
  }