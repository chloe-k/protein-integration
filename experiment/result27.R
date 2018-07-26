################################## Result 27 in GMR ############################################################
registerDoParallel(cores = 4)

# make RData after DRW

id_list <- c("27_0.2", "27_0.4", "27_0.6", "27_0.8", "27_0.9")
Gamma_list <- c(0.2, 0.4, 0.6, 0.8, 0.9)


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_27 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], prob = 0.001, Gamma = Gamma_list[i])
}

################################## Result 27_1 in GMR ############################################################
# Gamma = 0.2

id_list <- c("27_0.2_G", "27_0.2_M", "27_0.2_R", "27_0.2_GM", "27_0.2_GP", "27_0.2_MP", "27_0.2_GMP")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


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

