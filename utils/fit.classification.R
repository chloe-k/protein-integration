fit.classification <- function(y, samples, datapath, respath, profile_name, method = "DRW", pranking = "t-test", classifier = "glm", nFolds = 5, numTops=50, iter = 1){
  
  if(method == 'DRW'){
    
    load(file.path(datapath, paste(c("pathway_rank", profile_name, method, "RData"), collapse = '.')))

  } else if(method == 'gf'){
    
    load(file.path(datapath, paste(c("gene_rank", profile_name, method, "RData"), collapse = '.')))
    
  }
  
  
  
  
  Y <- rep(0,length(samples))
  Y[y[[1]]] <- 1
  Y=as.factor(Y)
  
  rankn_feats <- names(stats_feats)[1:numTops]
  result <- train(X[,rankn_feats], Y, trControl=trainControl(method="repeatedcv", number=nFolds, repeats = iter, returnResamp = "all"), method=classifier, family=binomial())
  
  return(result)
  
  
}