################################## Result 27 in GMR ############################################################

registerDoParallel(cores = 4)

# make RData after DRW

id_list <- c("27_0.2", "27_0.4", "27_0.6", "27_0.8", "27_0.9")
Gamma_list <- c(0.2, 0.4, 0.6, 0.8, 0.9)


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

id_list <- c("27_1_0.2", "27_2_0.2", "27_3_0.2", "27_4_0.2", "27_5_0.2", "27_6_0.2", "27_7_0.2")
res_gmr_27_0.2 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.2)
}

id_list <- c("27_1_0.4", "27_2_0.4", "27_3_0.4", "27_4_0.4", "27_5_0.4", "27_6_0.4", "27_7_0.4")
res_gmr_27_0.4 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.4)
}

id_list <- c("27_1_0.6", "27_2_0.6", "27_3_0.6", "27_4_0.6", "27_5_0.6", "27_6_0.6", "27_7_0.6")
res_gmr_27_0.6 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.6)
}

id_list <- c("27_1_0.8", "27_2_0.8", "27_3_0.8", "27_4_0.8", "27_5_0.8", "27_6_0.8", "27_7_0.8")
res_gmr_27_0.8 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.8)
}

id_list <- c("27_1_0.9", "27_2_0.9", "27_3_0.9", "27_4_0.9", "27_5_0.9", "27_6_0.9", "27_7_0.9")
res_gmr_27_0.9 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = type_list[i], prob = 0.001, Gamma = 0.9)
}


for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_27_', i, '_LOOCV.RData'), collapse = ''))
}


################################## Result 27_1 in GMR ############################################################
# Gamma = 0.2


res_gmr <- list(res_pa_GMR_27_1_LOOCV, res_pa_GMR_27_2_LOOCV, res_pa_GMR_27_3_LOOCV, res_pa_GMR_27_4_LOOCV, 
                res_pa_GMR_27_5_LOOCV, res_pa_GMR_27_6_LOOCV, res_pa_GMR_27_7_LOOCV)

for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR'), collapse = '')
  write.SigFeatures(res_fit=res_gmr[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

# Plot for GMR model
title <- c("Result 27_GMR")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr, perf_min = perf_min-0.02, perf_max = perf_max+0.02)




################################## Result 27 in GMR_d ############################################################
# make RData after DRW

id_list <- c("27_p0.2g0.2", "27_p0.2g0.4", "27_p0.2g0.6", "27_p0.2g0.8",
             "27_p0.4g0.2", "27_p0.4g0.4", "27_p0.4g0.6", "27_p0.4g0.8",
             "27_p0.6g0.2", "27_p0.6g0.4", "27_p0.6g0.6", "27_p0.6g0.8",
             "27_p0.8g0.2", "27_p0.8g0.4", "27_p0.8g0.6", "27_p0.8g0.8")

prob_list <- c(0.2, 0.2, 0.2, 0.2,
               0.4, 0.4, 0.4, 0.4,
               0.6, 0.6, 0.6, 0.6,
               0.8, 0.8, 0.8, 0.8)

Gamma_list <- c(0.2, 0.4, 0.6, 0.8,
                0.2, 0.4, 0.6, 0.8,
                0.2, 0.4, 0.6, 0.8,
                0.2, 0.4, 0.6, 0.8)

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_d_27 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], prob = prob_list[i], Gamma = Gamma_list[i])
}


################################## Result 27_1 in GMR_d ############################################################
# prob = 0.2
# Gamma = 0.2

id_list <- c("27_1", "27_2", "27_3", "27_4", "27_5", "27_6", "27_7")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_1 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.2, Gamma = 0.2)
}


################################## Result 27_2 in GMR_d ############################################################
# prob = 0.2
# Gamma = 0.4

id_list <- c("27_2_G", "27_2_M", "27_2_R", "27_2_GM", "27_2_GR", "27_2_MR", "27_2_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_2 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.2, Gamma = 0.4)
}


################################## Result 27_3 in GMR_d ############################################################
# prob = 0.2
# Gamma = 0.6

id_list <- c("27_3_G", "27_3_M", "27_3_R", "27_3_GM", "27_3_GR", "27_3_MR", "27_3_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_3 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.2, Gamma = 0.6)
}

################################## Result 27_4 in GMR_d ############################################################
# prob = 0.2
# Gamma = 0.8

id_list <- c("27_4_G", "27_4_M", "27_4_R", "27_4_GM", "27_4_GR", "27_4_MR", "27_4_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_4 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.2, Gamma = 0.8)
}


################################## Result 27_5 in GMR_d ############################################################
# prob = 0.4
# Gamma = 0.2

id_list <- c("27_5_G", "27_5_M", "27_5_R", "27_5_GM", "27_5_GR", "27_5_MR", "27_5_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_5 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.4, Gamma = 0.2)
}



################################## Result 27_6 in GMR_d ############################################################
# prob = 0.4
# Gamma = 0.4

id_list <- c("27_6_G", "27_6_M", "27_6_R", "27_6_GM", "27_6_GR", "27_6_MR", "27_6_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_6 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.4, Gamma = 0.4)
}



################################## Result 27_7 in GMR_d ############################################################
# prob = 0.4
# Gamma = 0.6

id_list <- c("27_7_G", "27_7_M", "27_7_R", "27_7_GM", "27_7_GR", "27_7_MR", "27_7_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_7 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.4, Gamma = 0.6)
}



################################## Result 27_8 in GMR_d ############################################################
# prob = 0.4
# Gamma = 0.8

id_list <- c("27_8_G", "27_8_M", "27_8_R", "27_8_GM", "27_8_GR", "27_8_MR", "27_8_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_8 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.4, Gamma = 0.8)
}



################################## Result 27_9 in GMR_d ############################################################
# prob = 0.6
# Gamma = 0.2

id_list <- c("27_9_G", "27_9_M", "27_9_R", "27_9_GM", "27_9_GR", "27_9_MR", "27_9_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_9 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.6, Gamma = 0.2)
}



################################## Result 27_10 in GMR_d ############################################################
# prob = 0.6
# Gamma = 0.4

id_list <- c("27_10_G", "27_10_M", "27_10_R", "27_10_GM", "27_10_GR", "27_10_MR", "27_10_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_10 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.6, Gamma = 0.4)
}



################################## Result 27_11 in GMR_d ############################################################
# prob = 0.6
# Gamma = 0.6

id_list <- c("27_11_G", "27_11_M", "27_11_R", "27_11_GM", "27_11_GR", "27_11_MR", "27_11_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_11 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.6, Gamma = 0.6)
}



################################## Result 27_12 in GMR_d ############################################################
# prob = 0.6
# Gamma = 0.8

id_list <- c("27_12_G", "27_12_M", "27_12_R", "27_12_GM", "27_12_GR", "27_12_MR", "27_12_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_12 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.6, Gamma = 0.8)
}



################################## Result 27_13 in GMR_d ############################################################
# prob = 0.8
# Gamma = 0.2

id_list <- c("27_13_G", "27_13_M", "27_13_R", "27_13_GM", "27_13_GR", "27_13_MR", "27_13_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_13 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.8, Gamma = 0.2)
}



################################## Result 27_14 in GMR_d ############################################################
# prob = 0.8
# Gamma = 0.4

id_list <- c("27_14_G", "27_14_M", "27_14_R", "27_14_GM", "27_14_GR", "27_14_MR", "27_14_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_14 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.8, Gamma = 0.4)
}



################################## Result 27_15 in GMR_d ############################################################
# prob = 0.8
# Gamma = 0.6

id_list <- c("27_15_G", "27_15_M", "27_15_R", "27_15_GM", "27_15_GR", "27_15_MR", "27_15_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_15 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.8, Gamma = 0.6)
}



################################## Result 27_16 in GMR_d ############################################################
# prob = 0.8
# Gamma = 0.8

id_list <- c("27_16_G", "27_16_M", "27_16_R", "27_16_GM", "27_16_GR", "27_16_MR", "27_16_GMR")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_27_16 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = type_list[i], prob = 0.8, Gamma = 0.8)
}

#############################################################################################
# write sigFeatures
for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '.RData'), collapse = ''))
}

res_gmr_d <- list(res_pa_GMR_d_27_2_G, res_pa_GMR_d_27_2_M, res_pa_GMR_d_27_2_R,
                  res_pa_GMR_d_27_2_GM, res_pa_GMR_d_27_2_GR, res_pa_GMR_d_27_2_MR,
                  res_pa_GMR_d_27_2_GMR)

profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Entrez)")
for(i in 1:length(id_list)){
  result_name <- paste(c('result',id_list[i],'_GMR_d'), collapse = '')
  write.SigFeatures(res_fit=res_gmr_d[[i]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)
}

##########################################################################################
########################################################################
# Plot for GMR_d (G)
title <- c("Result 27 GMR_d(G) Heatmap")
id_list <- c("27_1", "27_2_G", "27_3_G", "27_4_G",
             "27_5_G", "27_6_G", "27_7_G", "27_8_G", 
             "27_9_G", "27_10_G", "27_11_G", "27_12_G",
             "27_13_G", "27_14_G", "27_15_G", "27_16_G")

prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), each = 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)

res_gmr_d_27_G <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '.RData'), collapse = '')))
  res_gmr_d_27_G <- c(res_gmr_d_27_G, list(model))
}

perf_heatmap(title, res_gmr_d_27_G, prob_list = prob_list, Gamma_list = Gamma_list)

########################################################################
# Plot for GMR_d (M)
title <- c("Result 27 GMR_d(M) Heatmap")
id_list <- c("27_2", "27_2_M", "27_3_M", "27_4_M",
             "27_5_M", "27_6_M", "27_7_M", "27_8_M", 
             "27_9_M", "27_10_M", "27_11_M", "27_12_M",
             "27_13_M", "27_14_M", "27_15_M", "27_16_M")

prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), each = 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)

res_gmr_d_27_M <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '.RData'), collapse = '')))
  res_gmr_d_27_M <- c(res_gmr_d_27_M, list(model))
}

perf_heatmap(title, res_gmr_d_27_M, prob_list = prob_list, Gamma_list = Gamma_list)


########################################################################
# Plot for GMR_d (R)
title <- c("Result 27 GMR_d(R) Heatmap")
id_list <- c("27_3", "27_2_R", "27_3_R", "27_4_R",
             "27_5_R", "27_6_R", "27_7_R", "27_8_R", 
             "27_9_R", "27_10_R", "27_11_R", "27_12_R",
             "27_13_R", "27_14_R", "27_15_R", "27_16_R")

prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), each = 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)

res_gmr_d_27_R <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '.RData'), collapse = '')))
  res_gmr_d_27_R <- c(res_gmr_d_27_R, list(model))
}

perf_heatmap(title, res_gmr_d_27_R, prob_list = prob_list, Gamma_list = Gamma_list)


########################################################################
# Plot for GMR_d (GM)
title <- c("Result 27 GMR_d(GM) Heatmap")
id_list <- c("27_4", "27_2_GM", "27_3_GM", "27_4_GM",
             "27_5_GM", "27_6_GM", "27_7_GM", "27_8_GM", 
             "27_9_GM", "27_10_GM", "27_11_GM", "27_12_GM",
             "27_13_GM", "27_14_GM", "27_15_GM", "27_16_GM")

prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), each = 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)

res_gmr_d_27_GM <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '.RData'), collapse = '')))
  res_gmr_d_27_GM <- c(res_gmr_d_27_GM, list(model))
}

perf_heatmap(title, res_gmr_d_27_GM, prob_list = prob_list, Gamma_list = Gamma_list)


########################################################################
# Plot for GMR_d (GR)
title <- c("Result 27 GMR_d(GR) Heatmap")
id_list <- c("27_5", "27_2_GR", "27_3_GR", "27_4_GR",
             "27_5_GR", "27_6_GR", "27_7_GR", "27_8_GR", 
             "27_9_GR", "27_10_GR", "27_11_GR", "27_12_GR",
             "27_13_GR", "27_14_GR", "27_15_GR", "27_16_GR")

prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), each = 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)

res_gmr_d_27_GR <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '.RData'), collapse = '')))
  res_gmr_d_27_GR <- c(res_gmr_d_27_GR, list(model))
}

perf_heatmap(title, res_gmr_d_27_GR, prob_list = prob_list, Gamma_list = Gamma_list)


########################################################################
# Plot for GMR_d (MR)
title <- c("Result 27 GMR_d(MR) Heatmap")
id_list <- c("27_6", "27_2_MR", "27_3_MR", "27_4_MR",
             "27_5_MR", "27_6_MR", "27_7_MR", "27_8_MR", 
             "27_9_MR", "27_10_MR", "27_11_MR", "27_12_MR",
             "27_13_MR", "27_14_MR", "27_15_MR", "27_16_MR")

prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), each = 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)

res_gmr_d_27_MR <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '.RData'), collapse = '')))
  res_gmr_d_27_MR <- c(res_gmr_d_27_MR, list(model))
}

perf_heatmap(title, res_gmr_d_27_MR, prob_list = prob_list, Gamma_list = Gamma_list)


########################################################################
# Plot for GMR_d (GMR)
title <- c("Result 27 GMR_d(GMR) Heatmap")
id_list <- c("27_7", "27_2_GMR", "27_3_GMR", "27_4_GMR",
             "27_5_GMR", "27_6_GMR", "27_7_GMR", "27_8_GMR", 
             "27_9_GMR", "27_10_GMR", "27_11_GMR", "27_12_GMR",
             "27_13_GMR", "27_14_GMR", "27_15_GMR", "27_16_GMR")

prob_list <- rep(c(0.2, 0.4, 0.6, 0.8), each = 4)
Gamma_list <- rep(c(0.2, 0.4, 0.6, 0.8), 4)

res_gmr_d_27_GMR <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_', id_list[i], '.RData'), collapse = '')))
  res_gmr_d_27_GMR <- c(res_gmr_d_27_GMR, list(model))
}

perf_heatmap(title, res_gmr_d_27_GMR, prob_list = prob_list, Gamma_list = Gamma_list)


##################################GMR#########################################
res_gmr_27_0.2 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], prob = 0.001, Gamma = 0.2, type_used = type_list[i])
}

for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_', id_list[i], '_LOOCV.RData'), collapse = ''))
}



#GMR
title <- c("Result 27 GMR(Gamma = 0.2) ")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr_27_0.2, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_27_0.2, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_gmr_27_0.2, perf_min = perf_min-0.01, perf_max = perf_max+0.01)
perf_lineplot(fname_res = 'result/res_loocv_tuneK.txt', perf_min=55, perf_max=95)

################################## Result 27_2 in GMR ############################################################
# Gamma = 0.4

id_list <- c("27_0.4_G", "27_0.4_M", "27_0.4_R", "27_0.4_GM", "27_0.4_GP", "27_0.4_MP", "27_0.4_GMP")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_27_0.4 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], prob = 0.001, Gamma = 0.4, type_used = type_list[i])
}

for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_', id_list[i], '.RData'), collapse = ''))
}

res_gmr_27_0.4 <- list(res_pa_GMR_27_0.4_G, res_pa_GMR_27_0.4_M, res_pa_GMR_27_0.4_R,
                       res_pa_GMR_27_0.4_GM, res_pa_GMR_27_0.4_GP, res_pa_GMR_27_0.4_MP,
                       res_pa_GMR_27_0.4_GMP)

#GMR
title <- c("Result 27 GMR(Gamma = 0.4) ")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr_27_0.4, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_27_0.4, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr_27_0.4, perf_min = perf_min-0.01, perf_max = perf_max+0.01)

################################## Result 27_3 in GMR ############################################################
# Gamma = 0.6

id_list <- c("27_0.6_G", "27_0.6_M", "27_0.6_R", "27_0.6_GM", "27_0.6_GP", "27_0.6_MP", "27_0.6_GMP")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_27_0.6 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], prob = 0.001, Gamma = 0.6, type_used = type_list[i])
}

#GMR
title <- c("Result 27 GMR(Gamma = 0.6) ")
xlabs <- c("G", "M", "R", "GM", "GR", "MR", "GMR")

perf_min <- min(sapply(X = res_gmr_27_0.6, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_27_0.6, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr_27_0.6, perf_min = perf_min-0.01, perf_max = perf_max+0.01)

################################## Result 27_4 in GMR ############################################################
# Gamma = 0.8

id_list <- c("27_0.8_G", "27_0.8_M", "27_0.8_R", "27_0.8_GM", "27_0.8_GP", "27_0.8_MP", "27_0.8_GMP")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_27_0.8 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], prob = 0.001, Gamma = 0.8, type_used = type_list[i])
}


################################## Result 27_5 in GMR ############################################################
# Gamma = 0.9

id_list <- c("27_0.9_G", "27_0.9_M", "27_0.9_R", "27_0.9_GM", "27_0.9_GP", "27_0.9_MP", "27_0.9_GMP")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_27_0.9 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], prob = 0.001, Gamma = 0.9, type_used = type_list[i])
}

