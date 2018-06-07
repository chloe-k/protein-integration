unify_geneSymbol <- function(datapath){
  
  gdacpath <- file.path(datapath, 'BRCA_GDAC')
  
  # read gene symbol mapping table
  sym_rna <- read.delim('data/EDA/RNAseq(checker).txt', header = T, stringsAsFactors = F)
  sym_met <- read.delim('data/EDA/methyl(checker).txt', header = T, stringsAsFactors = F)
  sym_rppa <- read.delim('data/EDA/rppa(checker).txt', header = T, stringsAsFactors = F)
  sym_ppi <- read.delim('data/EDA/ppi(checker).txt', header = T, stringsAsFactors = F)
  sym_pathG <- read.delim('data/EDA/directGraph(checker).txt', header = T, stringsAsFactors = F)
  
  # merge gene symbol mapping table
  sym_map <- rbind(sym_rna, sym_met, sym_rppa, sym_ppi, sym_pathG)
  sym_map <- sym_map[!duplicated(sym_map),]
  
  # remove Approved name, HGNC ID, Location column
  # now sym_map has only Input, Match type, Approved symbol column
  sym_map <- sym_map[,-(4:6)]
  
  # remove rows which have Match type == Approved symbol
  sym_map <- sym_map[-which(sym_map$Match.type == 'Approved.symbol'),]
  
  # synonym processsssss!!!
  
  
  # RNAseq
  rna_sample <- names(read.csv(file.path(gdacpath, 'brca_rnaseq_qc.csv'), nrows = 1))
  rnaseq <- read.csv(file.path(gdacpath, 'brca_rnaseq_qc.csv'), header=F, row.names=1, col.names = rna_sample, skip=17)
  rnaseq$gene_id <- c()
  
  hgnc_rna <- mapping_rule(rownames(rnaseq), sym_map)
  
  # Methyl
  methyl <- read.csv(file.path(gdacpath, 'brca_methylation.csv'), header=T, row.names=1)
  
  # RPPA
  rppa <- read.delim(file.path(gdacpath, 'brca_rppa.txt'), header = T, stringsAsFactors = F, row.names = 1)
  
  # directed PPI
  ppi <- read.delim(file.path(datapath, 'dPPI_PathwayCommons9.sif'), header = T, stringsAsFactors = F, row.names = 1)
  
  # KEGG pathway graph
  #pathGraph <- 
  
}

mapping_rule <- function(gene, sym_map){
  
  for(i in 1:length(gene)){
    if(gene[i] %in% sym_map$Input){
      if(sym_map$Match.type == 'Previous symbol'){
        
      }else if(sym_map$Match.type == 'Synonyms'){
        
      }else if(sym_map$Match.type == 'Unmatched' || sym_map$Match.type == 'Entry withdrawn'){
        
      }
    }else{
      # In case that input gene symbol is not matched in HGNC gene symbol
      gene[i] <- 'removed'
    }
  }
  
 
  
  
  return(gene)
}

