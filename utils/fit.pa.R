fit.pa <- function(pathActivity, good, poor, profile_name, method, pranking, desc, 
                   numTops=50, iter=1, classifier = "glm", nFolds = 5,respath, DApath=NULL) {
  
  desc <- c(profile_name, method, "txt")
  
  # rank pathway activities
  # ranking = t-test / DA
  fname_rank = file.path(respath, paste(c("pathway_rank", pranking, desc), collapse = '.'))
  
  y = list(match(good, colnames(pathActivity)), match(poor, colnames(pathActivity)))
  
  stats_feats <- rankPathActivity(pathActivity = pathActivity, y = y, 
                                  ranking = pranking, fname=fname_rank, DApath=DApath)
  
  # write pathway ranking
  write.table(x=matrix(stats_feats, nrow=length(stats_feats), dimnames=list(names(stats_feats),"rank")),
              file=fname_rank, sep="\t", row.names=T, col.names=T)
  
  print('rank pathway activities done..')
  
  # perform 5-fold cross validation on logistic regression model
  library(caret)
  
  X <- t(pathActivity)
  Y <- rep(0,length(colnames(pathActivity)))
  Y[y[[1]]] <- 1
  Y=as.factor(Y)
  
  # classifier : glm / svmLinear
  acc <- c()
  for(k in seq(5,100,by=5)) {
    rankn_feats <- names(stats_feats)[1:k]
    set.seed(111)
    model <- train(X[,rankn_feats], Y, trControl=trainControl(method="repeatedcv", number=nFolds, repeats = iter), method=classifier, family=binomial())
    # model <- train(X[,rankn_feats], Y, method=classifier, trControl = trainControl(method="LOOCV", repeats = iter), family=binomial())
    
    acc <- c(acc, model$results$Accuracy)
  }
  
  df <- data.frame(k=seq(5,100,by=5), accuracy=acc)
  write.table(x=df,file = file.path(respath, paste(c("res_accuracy_tuneK", desc), collapse = '.')), row.names = F,quote = F)
  
  rankn_feats <- names(stats_feats)[1:df$k[which.max(df$accuracy)]]
  set.seed(111)
  return(train(X[,rankn_feats], Y, trControl=trainControl(method="repeatedcv", number=nFolds, repeats = iter), method=classifier, family=binomial()))
  # return(train(X[,rankn_feats], Y, method=classifier, trControl = trainControl(method="LOOCV", repeats = iter), family=binomial()))
}