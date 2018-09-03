fit.classification <- function(y, samples, id, datapath, respath, profile_name, method = "DRW", pranking = "t-test", classifier = "rf", nFolds = 5, numTops=50, iter = 10){
  
  desc <- c(id, method, "txt")
  
  if(method == 'DRW'){
    
    # iDRW
    load(file.path(datapath, paste(c("pathway_rank", id, profile_name, method, "RData"), collapse = '.')))
    
    # iDRW + DA
    # load(file.path(datapath, paste(c("pathway_rank", id, profile_name, method, pranking, "RData"), collapse = '.')))
    
  } else if(method == 'gf'){
    
    load(file.path(datapath, paste(c("gene_rank", id, profile_name, method, "RData"), collapse = '.')))
    
  }

  msg <- paste(c(id," classification start!"), collapse = '')
  print(msg)
  
  Y <- rep(0,length(samples))
  Y[y[[1]]] <- 1
  Y=as.factor(Y)
  
  ###### Setting for tuning hyperparameter and model evaluation #################
  set.seed(111)
  
  # dim(X)[2] - The number of pathway
  numTops = dim(X)[2]/2
  
  # Search for Top N pathway
  toppath <- file.path(respath, paste(c("res_accuracy_tuneK", desc), collapse = '.'))
  
  if(!file.exists(toppath)){
    acc <- c()
    trControl <- trainControl(method = "repeatedcv", number = nFolds, repeats = iter)
    
    for(k in seq(5,numTops,by=5)) {
      rankn_feats <- names(stats_feats)[1:k]
      set.seed(111)
      
      model <- train(x = X[,rankn_feats], y = Y, method=classifier, metric = "Accuracy", 
                     trControl=trControl, importance = TRUE)
      
      
      acc <- c(acc, max(model$results$Accuracy))
    }
    
    
    df <- data.frame(k=seq(5,numTops,by=5), accuracy=acc)
    write.table(x=df,file = toppath, row.names = F,quote = F)
  }
  else{
    df <- read.delim(file = toppath, header = T, sep = '')
  }
  
  print('Getting top N pathways is done..')
  
  # Model evaluation with top N pathway
  set.seed(111)
  rankn_feats <- names(stats_feats)[1:df$k[which.max(df$accuracy)]]
  nFolds = 5
  iter = 10
  
  # trControl <- trainControl(method = "repeatedcv", number = nFolds, repeats = iter)
  trControl <- trainControl(method = "LOOCV")
  
  # result <- train(X[,rankn_feats], Y, trControl=trainControl(method="repeatedcv", number=nFolds, repeats = iter, returnResamp = "all"), method=classifier, family=binomial())
  # result <- train(X[,rankn_feats], Y, trControl=trainControl(method="LOOCV"), method=classifier, family=binomial())
  result <- train(x = X[,rankn_feats], y = Y, method=classifier, metric = "Accuracy", 
                  trControl=trControl, importance = TRUE)
  
  print('classification is done')
  
  return(result)
  
  
}