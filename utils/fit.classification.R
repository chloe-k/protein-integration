fit.classification <- function(y, samples, id, datapath, respath, profile_name, method = "DRW", pranking = "t-test", classifier = "rf", nFolds = 5, numTops=50, iter = 10){
  
  if(method == 'DRW'){
    
    load(file.path(datapath, paste(c("pathway_rank", id, profile_name, method, "RData"), collapse = '.')))
    
  } else if(method == 'gf'){
    
    load(file.path(datapath, paste(c("gene_rank", id, profile_name, method, "RData"), collapse = '.')))
    
  }
  
  
  Y <- rep(0,length(samples))
  Y[y[[1]]] <- 1
  Y=as.factor(Y)
  
  ###### Setting for tuning hyperparameter and model evaluation #################
  set.seed(111)
  
  # dim(X)[2] - The number of pathway
  numTops = 100
  if(dim(X)[2] < 100)
    numTops = dim(X)[2]
  
  # 5-Fold CV with 10 iters
  trControl <- trainControl(method = "cv", number = nFolds, search = "grid", repeats = iter)
  
  # LOOCV
  # trControl <- trainControl(method="LOOCV", repeats = iter)
  
  ##############################################################################
  
  # Search best mtry(Hyperparameter for Random Forest)
  tuneGrid <- expand.grid(.mtry = c(seq(5,50,by=5)))
  
  rf_mtry <- train(x = X, y = Y, method = "rf", metric = "Accuracy",
                   tuneGrid = tuneGrid, trControl = trControl, importance = TRUE)
  
  opt_mtry <- rf_mtry$bestTune$mtry
  
  
  # Search for Top N pathway
  acc <- c()
  tuneGrid <- expand.grid(.mtry = opt_mtry)
  
  for(k in seq(5,numTops,by=5)) {
    rankn_feats <- names(stats_feats)[1:k]
    set.seed(111)
    model <- train(x = X[,rankn_feats], y = Y, method=classifier, metric = "Accuracy", 
                   tuneGrid = tuneGrid, trControl=trControl, importance = TRUE)

    acc <- c(acc, model$results$Accuracy)
  }
  
  df <- data.frame(k=seq(5,numTops,by=5), accuracy=acc)
  write.table(x=df,file = file.path(respath, paste(c("res_accuracy_tuneK", desc), collapse = '.')), row.names = F,quote = F)
  
  
  # Model evaluation with top N pathway
  set.seed(111)
  rankn_feats <- names(stats_feats)[1:df$k[which.max(df$accuracy)]]
  trControl <- trainControl(method = "LOOCV")
  
  # result <- train(X[,rankn_feats], Y, trControl=trainControl(method="repeatedcv", number=nFolds, repeats = iter, returnResamp = "all"), method=classifier, family=binomial())
  # result <- train(X[,rankn_feats], Y, trControl=trainControl(method="LOOCV"), method=classifier, family=binomial())
  result <- train(x = X[,rankn_feats], y = Y, method=classifier, metric = "Accuracy", 
                  tuneGrid = tuneGrid, trControl=trControl, importance = TRUE)
  
  return(result)
  
  
}