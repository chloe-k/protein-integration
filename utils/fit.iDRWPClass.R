fit.iDRWPClass <-
  function(x, y, testStatistic, profile_name, globalGraph = NULL, datapath, respath, pathSet,
           method = "DRW", samples, pranking = "t-test", mode = "GMP",
           classifier = "glm", nFolds = 5, numTops=50, id="Result0", prob=0.8,
           iter = 1, Gamma=0.3, AntiCorr = FALSE, DEBUG=TRUE) {
    
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
      if(mode == "GMP"){
        
        # get W0 of G & M 
        gm <- globalGraph[[1]] %du% globalGraph[[2]]
        W0 <- getW0(list(gene_weight[[1]], gene_weight[[2]]), gm)
        
        # get W0 of P
        p_W0 <- diffus_ppi(datapath = datapath, gene_weight = gene_weight[[3]], ppi = globalGraph[[3]], prob = prob)
        
        # concatenate W0 and p_W0
        W0 <- c(W0, p_W0, use.names = TRUE)
        if(DEBUG) cat('Getting W0 is done...')
        
        gmp <- gm %du% globalGraph[[3]]
        # get adjacency matrix of the (integrated) gene-gene graph
        # wpath <- file.path(datapath, paste(c(mode,"W","RData"), collapse = '.'))
        # if(!file.exists(wpath)){
          W = getW(datapath = datapath, G = gmp, gene_weight = gene_weight, mode = mode)
        # }
        # W = get(load(wpath))
      } 
      else{
        W0 <- getW0(gene_weight, globalGraph)
        # W0 <- diffus_ppi(datapath = datapath, gene_weight = gene_weight, ppi = globalGraph, prob = prob)
        if(DEBUG) cat('Getting W0 is done...')
        
        # get adjacency matrix of the (integrated) gene-gene graph
        # wpath <- file.path(datapath, paste(c(mode,"W","RData"), collapse = '.'))
        # if(!file.exists(wpath)){
          W = getW(datapath = datapath, G = globalGraph, gene_weight = gene_weight, mode = mode)
        # }
        # W = get(load(wpath))
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
      
      # if(!file.exists(pApath)){
        pA <- getPathActivity(x = x, pathSet = pathSet, w = vertexWeight, vertexZP = x_stats, 
                              method = method, fname = fname_profile, rows = samples)
        
        save(pA, file=pApath)
      # }else{
        # pA <- get(load(file = file.path(pApath)))
      # }

      
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
    #-----------------Classification-----------------#
    # perform 5-fold cross validation on logistic regression model
    
    # Y <- rep(0,length(samples))
    # Y[y[[1]]] <- 1
    # Y=as.factor(Y)
    
    # classifier : glm / svmLinear
    
    # acc <- c()
    # for(k in seq(5,100,by=5)) {
    #   rankn_feats <- names(stats_feats)[1:k]
    #   set.seed(111)
    #   #model <- train(X[,rankn_feats], Y, trControl=trainControl(method="repeatedcv", number=nFolds, repeats = iter), method=classifier, family=binomial())
    #   model <- train(X[,rankn_feats], Y, method=classifier, trControl = trainControl(method="LOOCV", repeats = iter), family=binomial())
    #   acc <- c(acc, model$results$Accuracy)
    # }
    # 
    # df <- data.frame(k=seq(5,100,by=5), accuracy=acc)
    # write.table(x=df,file = file.path(respath, paste(c("res_accuracy_tuneK", desc), collapse = '.')), row.names = F,quote = F)
    # 
    # rankn_feats <- names(stats_feats)[1:df$k[which.max(df$accuracy)]]
    # set.seed(111)
    # rankn_feats <- names(stats_feats)[1:numTops]
    # return(train(X[,rankn_feats], Y, trControl=trainControl(method="repeatedcv", number=nFolds, repeats = iter, returnResamp = "all"), method=classifier, family=binomial()))
    
  }