
#10628*228
cptac <- read.csv(file="iDRW/data/CPTAC/TCGA_Breast_BI_Proteome.itraq.tsv", sep="\t", header = T, row.names = 1)
#82720*74
cbio <- read.csv(file="data/brca/tcga/data_protein_quantification_Zscores.txt",sep="\t", header = T, row.names = 1)
#17673*870
rna <- read.csv(file="GIW/iDRW/data/BRCA_GDAC/brca_rnaseq_qc.csv")
#17037*869
met <- read.csv(file="GIW/iDRW/data/BRCA_GDAC/brca_methylation.csv")
#PPI
ppi <- read.csv(file="GIW/PathwayCommons9.All.hgnc.sif", sep="\t", header=F, stringsAsFactors = F)


#remove row with all na 
#50675*74
cbio <- cbio[rowSums(is.na(cbio))!=ncol(cbio),]

#sample
#rna_sample : 868
rna_sample <- colnames(rna)
rna_sample <- rna_sample[-1:-2]
#met_sample : 868
met_sample <- colnames(met)
met_sample <- met_sample[-1]
#cptac_sample : 105
cptac_sample <- colnames(cptac)
cptac_sample <- cptac_sample[-1:-7]
cptac_sample <- cptac_sample[-length(cptac_sample)+6:-length(cptac_sample)]
cptac_sample <- substring(cptac_sample,1,10)
cptac_sample <- paste("TCGA.",cptac_sample)
cptac_sample <- gsub(" ","",cptac_sample)
cptac_sample <- unique(cptac_sample)
#cbio : 74
cbio_sample <- colnames(cbio)

#rna gene : 17657
rna_gene <- rna[,1]
rna_gene <- rna_gene[-1:-16]
#met gene : 17037
met_gene <- met[,1]
#cptac gene : 10625
cptac_gene <- rownames(cptac)
cptac_gene <- cptac_gene[-1:-3]
#cbio gene : 11113
cbio_gene <- rownames(cbio)
cbio_gene_tmp <- strsplit(cbio_gene,"\\|")
cbio_gene <- sapply(cbio_gene_tmp,"[",1)
cbio_gene <- unique(cbio_gene)
#ppi gene : 31693
ppi_gene <- c(ppi[,1])
ppi_gene <- c(ppi_gene,ppi[,3])
ppi_gene <- unique(ppi_gene)

#sample
#rna & methyl : 568
rm <- length(intersect(rna_sample,met_sample))

#rna & methyl & cptac : 40 
rmc <- length(intersect(intersect(rna_sample,met_sample),cptac_sample))

#rna & methyl & cbio : 28
rmcb <- length(intersect(intersect(rna_sample,met_sample),cbio_sample))

#gene
#rna & methyl : 16454
grm <- length(intersect(rna_gene,met_gene))

#rna & methyl & cptac : 9044
grmc <-length(intersect(rna_gene,intersect(met_gene,cptac_gene)))

#rna & methyl & cbio : 9668
grmcb <- length(intersect(rna_gene,intersect(met_gene,cbio_gene)))

#rna & methyl & ppi : 14360
grmp <- length(intersect(rna_gene,intersect(met_gene,ppi_gene)))

#TCGA-BRCA 
tcga <- read.csv(file="most-affected-cases-table.tsv", sep="\t", header = T)
tcga_sample <- tcga$Case.ID
cptac_sample_tmp <- gsub("\\.","-",cptac_sample)
cptac_sample_tmp <- substring(cptac_sample_tmp,1,12)


#dim(gene,sample)
#cptac : 9044*40
#cbio : 9668*28
#rppa : 184*376
