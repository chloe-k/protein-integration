################################## Result 18 in GM ############################################################
num_cores <- detectCores()/2
registerDoParallel(cores = 4)

id_list <- c("18_1", "18_2", "18_3", "18_4", "18_5",
             "18_6", "18_7", "18_8", "18_9", "18_10")

Gamma_list <- c(0, 0.2, 0.4, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gm_18 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GM_model(id=id_list[i], type_used = "gmp", prob = 0.001, Gamma = Gamma_list[i])
}

for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GM_18_', i, '_LOOCV.RData'), collapse = ''))
}

# Plot for GM models
res_gm <- list(res_pa_GM_18_1_LOOCV, res_pa_GM_18_2_LOOCV, res_pa_GM_18_3_LOOCV, res_pa_GM_18_4_LOOCV, res_pa_GM_18_5_LOOCV,
               res_pa_GM_18_6_LOOCV, res_pa_GM_18_7_LOOCV, res_pa_GM_18_8_LOOCV, res_pa_GM_18_9_LOOCV, res_pa_GM_18_10_LOOCV)

title <- c("Result 18_GM")
xlabs <- c("[g=0]", "[g=0.2]", "[g=0.4]", "[g=0.6]", "[g=0.7]", "[g=0.75]", "[g=0.8]", "[g=0.85]", "[g=0.9]", "[g=0.95]")

perf_min <- min(sapply(X = res_gm, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gm, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gm, perf_min = perf_min-0.01, perf_max = perf_max+0.01, NULL)




################################## Result 18 in GMR ############################################################
num_cores <- detectCores()/2
registerDoParallel(cores = num_cores)

id_list <- c("18_1", "18_2", "18_3", "18_4", "18_5",
             "18_6", "18_7", "18_8")

Gamma_list <- c(0, 0.2, 0.4, 0.6, 0.8,
                0.85, 0.9, 0.95)

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_18 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], type_used = "gmp", prob = 0.001, Gamma = Gamma_list[i])
}

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_18_', i, '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}

############################################## plot #######################################

# Plot for GMR models

res_gmr <- list(res_pa_GMR_18_1_LOOCV, res_pa_GMR_18_2_LOOCV, res_pa_GMR_18_3_LOOCV, res_pa_GMR_18_4_LOOCV, res_pa_GMR_18_5_LOOCV,
                res_pa_GMR_18_6_LOOCV, res_pa_GMR_18_7_LOOCV, res_pa_GMR_18_8_LOOCV)

title <- c("Result 18_GMR")
xlabs <- c("[g=0]", "[g=0.2]", "[g=0.4]", "[g=0.6]", "[g=0.8]", "[g=0.85]", "[g=0.9]", "[g=0.95]")
perf_min <- min(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_gmr, perf_min = perf_min-0.15, perf_max = perf_max+0.15)


################################## Result 18 in GMR_d ############################################################
num_cores <- 4
registerDoParallel(cores = num_cores)

id_list <- c("18_1", "18_2", "18_3", "18_4", "18_5",
             "18_6", "18_7", "18_8", "18_9", "18_10",
             "18_11", "18_12", "18_13", "18_14", "18_15",
             "18_16", "18_17", "18_18", "18_19", "18_20",
             "18_21", "18_22", "18_23", "18_24", "18_25",
             "18_26", "18_27", "18_28", "18_29", "18_30")


prob_list <- c(0.001, 0.001, 0.001, 0.001, 0.001,
               0.01, 0.01, 0.01, 0.01, 0.01,
               0.2, 0.2, 0.2, 0.2, 0.2,
               0.4, 0.4, 0.4, 0.4, 0.4,
               0.6, 0.6, 0.6, 0.6, 0.6,
               0.8, 0.8, 0.8, 0.8, 0.8)

Gamma_list <- c(0, 0.2, 0.4, 0.6, 0.8,
                0, 0.2, 0.4, 0.6, 0.8,
                0, 0.2, 0.4, 0.6, 0.8,
                0, 0.2, 0.4, 0.6, 0.8,
                0, 0.2, 0.4, 0.6, 0.8,
                0, 0.2, 0.4, 0.6, 0.8)


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")

res_gmr_d_18 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], type_used = "gmp", prob = prob_list[i], Gamma = Gamma_list[i])
}

for(i in 1:length(id_list)){
  load(paste(c('data/model/res_pa_GMR_d_18_', i, '_LOOCV.RData'), collapse = ''))
}

############################################## plot #######################################
# Plot for GMR_d models
res_gmr_d <- list(res_pa_GMR_d_18_1_LOOCV, res_pa_GMR_d_18_2_LOOCV, res_pa_GMR_d_18_3_LOOCV, res_pa_GMR_d_18_4_LOOCV, res_pa_GMR_d_18_5_LOOCV,
                  res_pa_GMR_d_18_6_LOOCV, res_pa_GMR_d_18_7_LOOCV, res_pa_GMR_d_18_8_LOOCV, res_pa_GMR_d_18_9_LOOCV, res_pa_GMR_d_18_10_LOOCV,
                  res_pa_GMR_d_18_11_LOOCV, res_pa_GMR_d_18_12_LOOCV, res_pa_GMR_d_18_13_LOOCV, res_pa_GMR_d_18_14_LOOCV, res_pa_GMR_d_18_15_LOOCV,
                  res_pa_GMR_d_18_16_LOOCV, res_pa_GMR_d_18_17_LOOCV, res_pa_GMR_d_18_18_LOOCV, res_pa_GMR_d_18_19_LOOCV, res_pa_GMR_d_18_20_LOOCV,
                  res_pa_GMR_d_18_21_LOOCV, res_pa_GMR_d_18_22_LOOCV, res_pa_GMR_d_18_23_LOOCV, res_pa_GMR_d_18_24_LOOCV, res_pa_GMR_d_18_25_LOOCV,
                  res_pa_GMR_d_18_26_LOOCV, res_pa_GMR_d_18_27_LOOCV, res_pa_GMR_d_18_28_LOOCV, res_pa_GMR_d_18_29_LOOCV, res_pa_GMR_d_18_30_LOOCV)

title <- c("Result 18_GMR_d")
xlabs <- c("[p=0.001,g=0]", "[p=0.001,g=0.2]", "[p=0.001,g=0.4]", "[p=0.001,g=0.6]", "[p=0.001,g=0.8]",
           "[p=0.01,g=0]", "[p=0.01,g=0.2]", "[p=0.01,g=0.4]", "[p=0.01,g=0.6]", "[p=0.01,g=0.8]",
           "[p=0.2,g=0]", "[p=0.2,g=0.2]", "[p=0.2,g=0.4]", "[p=0.2,g=0.6]", "[p=0.2,g=0.8]",
           "[p=0.4,g=0]", "[p=0.4,g=0.2]", "[p=0.4,g=0.4]", "[p=0.4,g=0.6]", "[p=0.4,g=0.8]",
           "[p=0.6,g=0]", "[p=0.6,g=0.2]", "[p=0.6,g=0.4]", "[p=0.6,g=0.6]", "[p=0.6,g=0.8]",
           "[p=0.8,g=0]", "[p=0.8,g=0.2]", "[p=0.8,g=0.4]", "[p=0.8,g=0.6]", "[p=0.8,g=0.8]")

perf_min <- min(sapply(X = res_gmr_d, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_gmr_d, FUN = function(x){max(x$results$Accuracy)}))
perf_facet_boxplot(title, xlabs, res_gmr_d, perf_min = perf_min-0.03, perf_max = perf_max+0.03, perf_max)

title <- c("Result 18 GMR_d Heatmap")
perf_heatmap(title, xlabs, res_gmr_d)