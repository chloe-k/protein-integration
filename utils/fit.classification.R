fit.classification <- function(y, samples, id, datapath, respath, profile_name, method = "DRW", pranking = "t-test", classifier = "rf", nFolds = 5, numTops=50, iter = 20){
  
  if(method == 'DRW'){
    
    load(file.path(datapath, paste(c("pathway_rank", id, profile_name, method, "RData"), collapse = '.')))

  } else if(method == 'gf'){
    
    load(file.path(datapath, paste(c("gene_rank", id, profile_name, method, "RData"), collapse = '.')))
    
  }
  
  
  
  
  Y <- rep(0,length(samples))
  Y[y[[1]]] <- 1
  Y=as.factor(Y)
  
  acc <- c()
  for(k in seq(5,dim(X)[2],by=5)) {
    rankn_feats <- names(stats_feats)[1:k]
    set.seed(111)
    model <- train(X[,rankn_feats], Y, trControl=trainControl(method="repeatedcv", number=nFolds, repeats = iter), method=classifier, family=binomial())
    # model <- train(X[,rankn_feats], Y, method=classifier, trControl = trainControl(method="LOOCV", repeats = iter), family=binomial())
    
    acc <- c(acc, model$results$Accuracy)
  }
  
  df <- data.frame(k=seq(5,dim(x)[2],by=5), accuracy=acc)
  write.table(x=df,file = file.path(respath, paste(c("res_accuracy_tuneK", desc), collapse = '.')), row.names = F,quote = F)
  
  rankn_feats <- names(stats_feats)[1:df$k[which.max(df$accuracy)]]
  set.seed(111)
  
  # result <- train(X[,rankn_feats], Y, trControl=trainControl(method="repeatedcv", number=nFolds, repeats = iter, returnResamp = "all"), method=classifier, family=binomial())
  result <- train(X[,rankn_feats], Y, trControl=trainControl(method="LOOCV"), method=classifier, family=binomial())
  
  return(result)
  
  
}