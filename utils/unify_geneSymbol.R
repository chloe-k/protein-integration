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
  sym_map <- sym_map[-which(sym_map$Match.type == 'Approved symbol'),]
  
  # remove rows which have Match type == Synonyms
  sym_map <- sym_map[-which(sym_map$Match.type == 'Synonyms'),]
  
  #####################################################################################################################
  
  # RNAseq
  rna_sample <- names(read.csv(file.path(gdacpath, 'brca_rnaseq_qc.csv'), nrows = 1))
  rnaseq <- read.csv(file.path(gdacpath, 'brca_rnaseq_qc.csv'), header=F, row.names=1, col.names = rna_sample, skip=17)
  rnaseq$gene_id <- c()
  
  hgnc_rna <- mapping_rule(rnaseq, sym_map)
  write.csv(hgnc_rna, file="HGNC_brca_rnaseq_qc.csv", sep=",")
  
  # Methyl
  methyl <- read.csv(file.path(gdacpath, 'brca_methylation.csv'), header=T, row.names=1)
  
  hgnc_met <- mapping_rule(methyl, sym_map)
  write.csv(hgnc_met, file="HGNC_brca_methylation.csv", sep=",")
  
  # RPPA
  rppa <- read.delim(file.path(gdacpath, 'brca_rppa.txt'), header = T, stringsAsFactors = F, row.names = 1)
  
  hgnc_rppa <- mapping_rule(rppa, sym_map)
  write.table(hgnc_rppa, file="HGNC_brca_rppa.txt", sep="\t")
  
  # directed PPI
  ppi <- read.delim(file.path(datapath, 'dPPI_PathwayCommons9.sif'), header = T, stringsAsFactors = F, row.names = 1)
  
  hgnc_ppi <- mapping_rule()
  
  # KEGG pathway graph
  #pathGraph <- 
  
}

mapping_rule <- function(data, sym_map){
  # This function performs mapping gene symbol to HGNC gene symbol
  
  cat('Checking profile dimension before mapping\n')
  cat(dim(data))
  cat('\n')
  
  gene <- rownames(data)
  
  # make mapping vector by using originnal gene name list of data
  for(i in 1:length(gene)){
    if(gene[i] %in% sym_map$Input){
      geneInsym_map <- sym_map[which(sym_map$Input == gene[i]),]
      if(geneInsym_map$Match.type == 'Previous symbol'){
        gene[i] <- geneInsym_map$Approved.symbol
      }else if(geneInsym_map$Match.type == 'Unmatched' || geneInsym_map$Match.type == 'Entry withdrawn'){
        gene[i] <- 'removed'
      }
    }
  }
  
  # remove rows which are not contained in HGNC symbol
  data <- data[-(which(gene == 'removed')),]
  gene <- gene[-which(gene == 'removed')]
  
  rownames(data) <- gene
  
  cat('Dimension after mapping\n')
  cat(dim(data))
  cat('\n')
  
  return(data)
}

mapping_rule_graph <- function(data, sym_map){
  # This function performs 3 steps
  # 1. extract unique graph node(gene) 
  # 2. make mapping table that consists of Before(unique graph node(gene)) and After(HGNC gene symbol)
  # 3. apply mapping table to data(Interaction data)
  
  cat('Dimension before mapping\n')
  cat(dim(data))
  cat('\n')
  
  gene <- data[,1]
  gene <- rbind(gene, data[,2])
  
  # step 1. all unique genes in Graph
  gene <- unique(gene)
  
  # step 2. make mapping table for graph interaction
  gene <- cbind(gene, "")
  colnames(gene) <- c("Before", "After")
  
  for(i in 1:length(gene)){
    if(gene$Before[i] %in% sym_map$Input){
      geneInsym_map <- sym_map[which(sym_map$Input == gene$Before[i]),]
      if(geneInsym_map$Match.type == 'Previous symbol'){
        gene$After[i] <- geneInsym_map$Approved.symbol
      }else if(geneInsym_map$Match.type == 'Unmatched' || geneInsym_map$Match.type == 'Entry withdrawn'){
        gene$After[i] <- 'removed'
      }
    }
  }
  
  # step 3. map HGNC gene symbol to gene symbol in graph

  data <- data.frame(lapply(X = data, 
                            FUN = function(x){
                              if(x %in% gene$Before){
                                gsub(x, gene$Before, gene$After)
                              }else x
                            }))

  cat('Dimension after mapping\n')
  cat(dim(data))
  cat('\n')
  
  return(data)
}

