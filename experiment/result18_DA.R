# DA was applied in Result 18
# 202.30.3.222 server

####################################### Result 18_DA ##############################################
method = "DRW"
AntiCorr = FALSE
###########################################GM#################################
id = "result18_3_GM"
profile_name <- c("rna(Entrez)", "meth(Entrez)")
da_weight_file <- "GM_compressed_data2.tsv"
pranking <- "DA"

pApath <- file.path(respath, paste(c("pA", id, profile_name, method, if(AntiCorr) "anticorr", "RData"), collapse = '.'))

load(file = pApath)

# rank pathway activities
# ranking = t-test / DA

fname_rank = file.path(datapath, paste(c("pathway_rank", pranking, profile_name, method, if(AntiCorr) "anticorr", "txt"), collapse = '.'))
DApath <- file.path(datapath, "DA_result", da_weight_file)
stats_feats <- rankPathActivity(ranking=pranking, fname=fname_rank, DApath=DApath)

X <- t(pA$pathActivity)

save(stats_feats, X, file=file.path(datapath, paste(c("pathway_rank", id, profile_name, method, pranking, "RData"), collapse = '.')))

# write pathway ranking
write.table(x=matrix(stats_feats, nrow=length(stats_feats), dimnames=list(names(stats_feats),"rank")),
            file=fname_rank, sep="\t", row.names=T, col.names=T)

res_pa_GM_18_DA <- fit.classification(y=y, samples = samples, id = id, datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "DA", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)

save(res_pa_GM_18_DA, file="data/model/res_pa_GM_18_DA.RData")


########################################### GMR #################################
id = "result18_3_GMR"
profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
da_weight_file <- "GMR_compressed_data2.tsv"
pranking <- "DA"

pApath <- file.path(respath, paste(c("pA", id, profile_name, method, if(AntiCorr) "anticorr", "RData"), collapse = '.'))

load(file = pApath)

# rank pathway activities
# ranking = t-test / DA

fname_rank = file.path(datapath, paste(c("pathway_rank", pranking, profile_name, method, if(AntiCorr) "anticorr", "txt"), collapse = '.'))
DApath <- file.path(datapath, "DA_result", da_weight_file)
stats_feats <- rankPathActivity(ranking=pranking, fname=fname_rank, DApath=DApath)

X <- t(pA$pathActivity)

save(stats_feats, X, file=file.path(datapath, paste(c("pathway_rank", id, profile_name, method, pranking, "RData"), collapse = '.')))

# write pathway ranking
write.table(x=matrix(stats_feats, nrow=length(stats_feats), dimnames=list(names(stats_feats),"rank")),
            file=fname_rank, sep="\t", row.names=T, col.names=T)

res_pa_GMR_18_DA <- fit.classification(y=y, samples = samples, id = id, datapath = datapath, respath = respath,
                                      profile_name = profile_name, method = "DRW", pranking = "DA", classifier = "rf",
                                      nFolds = 5, numTops=50, iter = 10)

save(res_pa_GMR_18_DA, file="data/model/res_pa_GMR_18_DA.RData")
