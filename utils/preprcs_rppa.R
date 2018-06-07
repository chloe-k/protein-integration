preprcs_rppa <- function(gdacpath){
  
  # non-unified gene was used
  #rawrppapath <- file.path(gdacpath, 'brca_rppa.txt')
  rawrppapath <- file.path(gdacpath, 'HGNC_brca_rppa.txt')
  
  # read rppa data (226*937)
  rppa <- read.delim(rawrppapath, header = T, stringsAsFactors = F, row.names = 1)
  
  multiGene <- vector()
  decomp_rppa <- rppa
  
  for(i in 1:length(rownames(decomp_rppa))){
    geneAntibody <- rownames(decomp_rppa)[i]
    gene <- strsplit(geneAntibody, "\\|")[[1]][1]
    geneSets <- strsplit(gene, " ")
    numGene <- length(geneSets[[1]])
    
    if(numGene > 1){
      multiGene <- append(multiGene, i)
      for(j in 1:numGene){
        decomp_rppa <- rbind(decomp_rppa, rep(decomp_rppa[i,], times=1))
        row.names(decomp_rppa)[length(decomp_rppa[,1])] <- paste(geneSets[[1]][j],"|",strsplit(geneAntibody,"\\|")[[1]][2],sep="")
      }
    }
  }
  decomp_rppa <- decomp_rppa[-multiGene,]
  
  # fill missing values with median(row) imputation
  library(mice)
  library(pracma)
  tail(md.pattern(decomp_rppa))
  
  imputed_rppa <- data.frame(apply(X=decomp_rppa,
                                   MARGIN = 1,
                                   FUN=function(x) {
                                     if(is.numeric(x))
                                       ifelse(is.na(x),median(x,na.rm=T),x)
                                     else x}))
  imputed_rppa <- as.data.frame(t(imputed_rppa))
  row.names(imputed_rppa) <- rownames(decomp_rppa)
  tail(md.pattern(imputed_rppa))
  
  #write.csv(imputed_rppa, file.path(gdacpath,'decomposed_rppa.csv'))
  write.csv(imputed_rppa, file.path(gdacpath,'HGNC_decomposed_rppa.csv'))
  
  # mean expression by gene
  # Because row name can not be duplicated, process mean aggregate in first column.
  rpparow <- sapply(strsplit(row.names(imputed_rppa), "\\|"), "[[", 1)
  imputed_rppa <- cbind(rpparow,imputed_rppa)
  
  mean_imputed_rppa <- aggregate(imputed_rppa[,-1], list(gene=imputed_rppa[,1]), FUN=mean)
  row.names(mean_imputed_rppa) <- mean_imputed_rppa[,1]
  mean_imputed_rppa <- mean_imputed_rppa[,-1]
  row.names(mean_imputed_rppa) <- trimws(rownames(mean_imputed_rppa), "b")
  
  
  #write.csv(mean_imputed_rppa, file.path(gdacpath,'mean_imputed_rppa.csv'))
  write.csv(mean_imputed_rppa, file.path(gdacpath,'HGNC_mean_imputed_rppa.csv'))
  
}