read_data <- function(year, datapath){
  
  gdacpath <- file.path(datapath, 'BRCA_GDAC')
  
  # read BRCA data
  methyl <- read.csv(file.path(gdacpath, 'brca_methylation.csv'), header=T, row.names=1)
  rnaseq <- read.csv(file.path(gdacpath, 'brca_rnaseq_qc.csv'), header=T, row.names=2)
  process_rppa(gdacpath)
  
  # read rppa data
  rppa <- read.csv(file.path(gdacpath, 'mean_imputed_rppa.csv'), header=T, row.names=1, stringsAsFactors = F)
  
  # process rppa
  row.names(rppa) <- substring(rownames(rppa),2) # 188 genes
  rppa_gene_id <- gene_name_id_map[rownames(rppa),]
  unmapped_rppa <- which(is.na(rppa_gene_id)) # 4 genes
  rppa <- rppa[-unmapped_rppa,]   # 184 genes
  
  # map gene name to id
  gene_name_id_map <- read.csv(file.path(gdacpath, 'gene_name_id_map'), skip=29,header=F, row.names=1)
  row.names(methyl) <- gene_name_id_map[rownames(methyl),]
  row.names(rppa) <- gene_name_id_map[rownames(rppa),]
  
  # remove gene name column (RNAseq)
  rnaseq$gene_name <- c()
  
  # differentiate RNAseq, methylation, RPPA genes
  row.names(rnaseq) <- paste("g", rownames(rnaseq), sep="")
  row.names(methyl) <- paste("m", rownames(methyl), sep="")
  row.names(rppa) <- paste("r", rownames(rppa), sep="")
  
  # read clinical information of BRCA patients
  clinical <- read.csv(file.path(datapath, 'brca_clinical.csv'), header=T, row.names=1)
  
  group_cut <- 365*year
  
  survival_all <- clinical$survival
  samples_all <- rownames(clinical)
  
  # extract overlapping samples
  samples_g <- colnames(rnaseq) # 868 samples
  samples_m <- colnames(methyl) # 868 samples
  samples_r <- colnames(rppa) # 937 samples
  samples <- intersect(samples_g, samples_m) # 568 samples
  
  samples_r <- gsub('\\_', '.', samples_r)
  samples_r <- substring(samples_r, 1,15)
  colnames(rppa) <- samples_r
  
  samples <- intersect(samples, samples_r) # 376 samples
  
  clinical <- clinical[samples,]
  
  # remove samples whose survival days were not recorded (NA) or wrongly so as negative values
  remove_samples <- which(is.na(clinical$survival))
  remove_samples <- c(remove_samples, which(clinical$survival<0))
  
  # remove samples whose survival days < 3 years and reported as "living"
  remove_samples <- c(remove_samples, which(clinical$survival<=group_cut & clinical$vital_status==1))
  
  samples <- samples[-remove_samples]
  clinical <- clinical[-remove_samples,]
  
  #total 376 samples
  rnaseq <- rnaseq[,samples]
  methyl <- methyl[,samples]
  rppa <- rppa[,samples]
  
  # split into 218 good / 247 poor group
  good_samples <- which(clinical$survival>group_cut)
  poor_samples <- which(clinical$survival<=group_cut)
  
  library(mice)
  tail(md.pattern(methyl)) # 5134 missing values
  
  # fill missing values with median imputation
  imputed_methyl <- data.frame(lapply(X=methyl,
                                      FUN=function(x) {
                                        if(is.numeric(x)) 
                                          ifelse(is.na(x),median(x,na.rm=T),x) 
                                        else x}))
  
  row.names(imputed_methyl) <- rownames(methyl)
  
  tail(md.pattern(imputed_methyl))
  
  
  # save as RData
  save(rnaseq, imputed_methyl, rppa, gene_name_id_map, clinical, samples, good_samples, poor_samples, file = file.path(datapath, 'data.RData'))
  
}