fit.classification <- function(y, samples, respath, profile_name, method = "DRW", pranking = "t-test", classifier = "glm", nFolds = 5, numTops=50, iter = 1){
  
  pA <- get(load(file.path(respath, paste(c("pA", profile_name, method, "RData"), collapse = '.'))))
  X <- t(pA$pathActivity)
  
  desc <- c(profile_name, method, "txt")
  fname_rank = file.path(respath, paste(c("pathway_rank", pranking, desc), collapse = '.'))
  stats_feats <- read.delim(file = fname_rank, header = T, sep = "\t", row.names = 1, stringsAsFactors = F)
  
  Y <- rep(0,length(samples))
  Y[y[[1]]] <- 1
  Y=as.factor(Y)
  
  rankn_feats <- rownames(stats_feats)[1:numTops]
  result <- train(X[,rankn_feats], Y, trControl=trainControl(method="repeatedcv", number=nFolds, repeats = iter, returnResamp = "all"), method=classifier, family=binomial())
  
  return(result)
  
  
}