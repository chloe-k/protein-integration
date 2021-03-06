count.sigFeatures <- function(res_fit, iter, nFolds) {
  # union significant pathway / gene features across 5 folds
  res <- res_fit[[4]]
  r <- res[[2]]
  p <- names(r)
  for (i in 2:iter*nFolds) {
    res <- res[[1]]
    r <- res[[2]]
    p <- union(p,names(r))
  }
  pathway_cnt <- data.frame(pathway_id=p, count=rep(0,length(p)), row.names = 1)
  
  # count significant pathway / gene features across 5 folds
  res <- res_fit[[4]]
  r <- res[[2]]
  pathway_cnt[names(r),] <- pathway_cnt[names(r),] + 1
  for (i in 2:iter*nFolds) {
    res <- res[[1]]
    r <- res[[2]]
    pathway_cnt[names(r),] <- pathway_cnt[names(r),] + 1
  }
  return(list(p, pathway_cnt))
}

# write significant pathway / gene features
write.SigFeatures <- function(res_fit, id, profile_name, method = "DRW", classifier = NULL, respath, AntiCorr = FALSE, da = FALSE) {
  
  pA <- get(load(file.path(respath, paste(c("pA", id, profile_name, method, if(AntiCorr) "anticorr", "RData"), collapse = '.'))))
  sigGeneset <- pA$sigGenes
  
  feats <- as.matrix(varImp(res_fit)$importance)
  feats <- feats[order(-feats[,1]),]
  
  genemap <- get.geneMapTable()
  
  #p <- substring(rownames(feats),2,6) 
  p <- rownames(feats) 
  # p <- names(feats) 
  
  pathway_name <- sapply(X = p, FUN = function(x) strsplit(keggGet(paste(c("hsa", x), collapse = ""))[[1]]$NAME, " - ")[[1]][1])
  
  desc <- c(id, profile_name, method, classifier, if(AntiCorr) "anticorr", if(da) "da","txt")
  fname_res <- file.path(respath, paste(c("sigPathway_genes", desc), collapse = '.'))
  
  sink(fname_res)
  for (i in 1:length(p)) {
    # sigGenes <- unlist(sigGeneset[p[i]],use.names = T)
    # cat(p[i], ';', as.character(pathway_name[p[i]]), ';', length(pathSet[[p[i]]]), ';', writeGeneSymbol(substring(sigGenes,2),substring(sigGenes,1,1)), '\n')
    sigGenes <- unlist(sigGeneset[p[i]],use.names = F)
    cat(p[i], ';', as.character(pathway_name[p[i]]), ';', length(pathSet[[p[i]]]), ';', writeGeneId(substring(sigGenes,2),genemap,substring(sigGenes,1,1)), '\n')
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

writeGeneId <- function(gene, genemap, gid){
  symbol <- as.character(genemap[gene,])
  for(i in 1:length(symbol)){
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

get.geneMapTable <- function() {
  genemap <- read.table('data/BRCA_GDAC/gene_name_id_map', sep = ',', row.names=2)
  return(genemap)
}