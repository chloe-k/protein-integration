read_data <- function(year, datapath){
  
  gdacpath <- file.path(datapath, 'BRCA_GDAC')
  
  # ------------------------------------- non-unified gene(HGNC gene symbol) was used --------------------------------------------#
  # read BRCA data
  methyl <- read.csv(file.path(gdacpath, 'brca_methylation.csv'), header=T, row.names=1)

  # remove unknown gene name(row)
  rna_sample <- names(read.csv(file.path(gdacpath, 'brca_rnaseq_qc.csv'), nrows = 1))
  rnaseq <- read.csv(file.path(gdacpath, 'brca_rnaseq_qc.csv'), header=F, row.names=1, col.names = rna_sample, skip=17)
  rnaseq$gene_id <- c()

  # imputation ppi node value by using RPPA data
  rppapath <- file.path(gdacpath, 'mean_imputed_rppa.csv')
  if(!file.exists(rppapath)) {
    print('preprocessed rppa profile does not exist, now preprocessing RPPA profile start')
    preprcs_rppa(gdacpath)
  }

  # read rppa data (188*937)
  rppa <- read.csv(file.path(gdacpath, 'mean_imputed_rppa.csv'), header=T, row.names=1, stringsAsFactors = F)
  # ------------------------------------- non-unified gene(HGNC gene symbol) was used --------------------------------------------#
  
  # ------------------------------------- unified gene(HGNC gene symbol) was used ------------------------------------------------#
  # # read BRCA data
  # methyl <- read.csv(file.path(gdacpath, 'HGNC_brca_methylation.csv'), header=T, row.names=1)
  # 
  # # remove unknown gene name(row)
  # rna_sample <- names(read.csv(file.path(gdacpath, 'HGNC_brca_rnaseq_qc.csv'), nrows = 1))
  # rnaseq <- read.csv(file.path(gdacpath, 'HGNC_brca_rnaseq_qc.csv'), header=F, row.names=1, col.names = rna_sample, skip=17)
  # rnaseq$gene_id <- c()
  # 
  # # imputation ppi node value by using RPPA data
  # rppapath <- file.path(gdacpath, 'HGNC_mean_imputed_rppa.csv')
  # if(!file.exists(rppapath)) {
  #   print('preprocessed rppa profile does not exist, now preprocessing RPPA profile start')
  #   preprcs_rppa(gdacpath)
  # }
  # 
  # # read rppa data (188*937)
  # rppa <- read.csv(file.path(gdacpath, 'HGNC_mean_imputed_rppa.csv'), header=T, row.names=1, stringsAsFactors = F)
  # ------------------------------------- unified gene(HGNC gene symbol) was used ------------------------------------------------#
  
  # differentiate RNAseq, methylation, RPPA genes
  row.names(rnaseq) <- paste("g", rownames(rnaseq), sep="")
  row.names(methyl) <- paste("m", rownames(methyl), sep="")
  row.names(rppa) <- paste("p", rownames(rppa), sep="")
  
  # read clinical information of BRCA patients
  clinical <- read.csv(file.path(gdacpath, 'brca_clinical.csv'), header=T, row.names=1)
  
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
  
  # split into 177 good / 199 poor group
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
  
  # non-unified gene symbol was used
  #save(rnaseq, imputed_methyl, rppa, clinical, samples, good_samples, poor_samples, file = file.path(datapath, 'data.RData'))
  
  # unified gene symbol(HGNC gene symbol) was used
  save(rnaseq, imputed_methyl, rppa, clinical, samples, good_samples, poor_samples, file = file.path(datapath, 'HGNC_unfy_data.RData'))
  
}