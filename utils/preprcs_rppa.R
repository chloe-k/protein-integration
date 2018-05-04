preprcs_rppa <- function(datapath, rawrppapath){
  rawrppapath <- file.path(datapath, 'brca_rppa.txt')
  
  # read rppa data (226*937)
  rppa <- read.delim(rawrppapath, header = T, stringsAsFactors = F, row.names = 1)
  
  #14
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
  tail(md.pattern(decomp_rppa))
  
  imputed_rppa <- data.frame(lapply(X=decomp_rppa,
                                    FUN=function(x) {
                                      if(is.numeric(x))
                                        ifelse(is.na(x),median(x,na.rm=T),x)
                                      else x}))
  
}