# write significant pathway / gene features
write.SigFeatures <- function(res_fit, profile_name, method = "DRW", classifier = NULL, respath, AntiCorr = FALSE, da = FALSE) {
  
  pA <- get(load(file.path(respath, paste(c("pA", profile_name, method, if(AntiCorr) "anticorr", "RData"), collapse = '.'))))
  sigGeneset <- pA$sigGenes
  
  feats <- as.matrix(varImp(res_fit)$importance)
  feats <- feats[order(-feats[,1]),]
  
  #p <- substring(names(feats),2,6)
  p <- names(feats)
  
  pathway_name <- sapply(X = p, FUN = function(x) strsplit(keggGet(paste(c("hsa", x), collapse = ""))[[1]]$NAME, " - ")[[1]][1])
  
  desc <- c(profile_name, method, classifier, if(AntiCorr) "anticorr", if(da) "da","txt")
  fname_res <- file.path(respath, paste(c("sigPathway_genes", desc), collapse = '.'))
  
  sink(fname_res)
  for (i in 1:length(p)) {
    sigGenes <- unlist(sigGeneset[p[i]],use.names = T)
    cat(p[i], ';', as.character(pathway_name[p[i]]), ';', length(pathSet[[p[i]]]), ';', writeGeneSymbol(substring(sigGenes,2),substring(sigGenes,1,1)), '\n')
  }
  sink.reset()
}



writeGeneSymbol <- function(symbol, gid) {
  for (i in 1:length(symbol)) {
    if(gid[i]=="g") symbol[i] <- paste(symbol[i], "(gene)",sep="")
    else if(gid[i]=="m") symbol[i] <- paste(symbol[i], "(meth)", sep="")
    else if(gid[i]=="p") symbol[i] <- paste(symbol[i], "(prot)", sep="")
  }
  return(symbol)
}

sink.reset <- function() {
  
  for (i in seq_len(sink.number())) {
    sink(NULL)
    sink(NULL, type = "message")
  }
  
}