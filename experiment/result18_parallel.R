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
perf_lineplot(title = title, xlabs = xlabs, res_models = res_models, perf_max = perf_max, perf_min = perf_min, Gamma_list = Gamma_list)

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

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_d_18_', i, '_LOOCV.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}

result_name <- paste(c('result',id_list[17],'_GMR_d'), collapse = '')
write.SigFeatures(res_fit=res_models[[17]], id = result_name, profile_name=profile_name, method="DRW", respath=respath)


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


############################################## plot #######################################
# Plot for comparison with GMR_26, GMR_d_26, GM_baseline
title <- c("Result 18")
g <- get(load('data/model/res_pa_G_11.RData'))
m <- get(load('data/model/res_pa_M_11.RData'))
r <- get(load('data/model/res_pa_R_11.RData'))
gm_rdc <- get(load('data/model/res_pa_GM_26.RData'))
gmr_rdc <- get(load('data/model/res_pa_GMR_26.RData'))
gmr_d_rdc <- get(load('data/model/res_pa_GMR_d_26.RData'))
gm <- get(load('data/model/res_pa_GM_18_3_LOOCV.RData'))
gmr <- get(load('data/model/res_pa_GMR_18_3_LOOCV.RData'))
gmr_d <- get(load('data/model/res_pa_GMR_d_18_17_LOOCV.RData'))
gr <- get(load('data/model/res_pa_GR_28_0.6_LOOCV.RData'))
gr_d <- get(load('data/model/res_pa_GR_d_29_16_LOOCV.RData'))

res_models <- list(g, m, r, gm_rdc, gmr_rdc, gmr_d_rdc, gm, gmr, gmr_d, gr, gr_d)

# write sigPathway (202.30.3.222)
# result_name <- paste(c('result18_3_GM'), collapse = '')
# write.SigFeatures(res_fit=gm, id = result_name, profile_name=profile_name, method="DRW", respath=respath)
# result_name <- paste(c('result18_3_GMR'), collapse = '')
# write.SigFeatures(res_fit=gmr, id = result_name, profile_name=profile_name, method="DRW", respath=respath)
# result_name <- paste(c('result18_3_GMR_d'), collapse = '')
# write.SigFeatures(res_fit=gmr_d, id = result_name, profile_name=profile_name, method="DRW", respath=respath)


xlabs <- c("DRW(G)", "DRW(M)", "DRW(P)",
           "iDRW(overlapped-GM)", "iDRW(overlapped-GMP)", "iDRW_pr(overlapped-GMP)",
           "iDRW(GM)", "iDRW(GMP)", "iDRW_pr(GMP)",
           "iDRW(GP)", "iDRW_pr(GP)")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.01, perf_max = perf_max+0.01)
perf_barplot(xlabs, res_models, perf_min-0.01, 0.75)



########## For single pathway profile ########
g <- get(load('data/model/res_pa_G_11.RData'))
m <- get(load('data/model/res_pa_M_11.RData'))
r <- get(load('data/model/res_pa_R_11.RData'))

res_models <- list(g, m, r)
xlabs <- c("DRW(G)", "DRW(M)", "DRW(P)")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_barplot(xlabs, res_models, 0.6, perf_max+0.01)


########## For reduced profile by extracting overlapped with RPPA protein (Result 26) ########
#380*610
gm_rdc <- get(load('data/model/res_pa_GM_26.RData'))
gmr_rdc <- get(load('data/model/res_pa_GMR_26.RData'))
gmr_d_rdc <- get(load('data/model/res_pa_GMR_d_26.RData'))

res_models <- list(gm_rdc, gmr_rdc, gmr_d_rdc)
# xlabs <- c(parse(text = paste('iDRW(G','^R','M','^R)', sep='')), parse(text = paste('iDRW(G','^R','M','^R','P)', sep='')), parse(text = paste('iDRW_pr(G','^R','M','^R','P)', sep='')))
xlabs <- c("iDRW(overlapped-GM)", "iDRW(overlapped-GMP)", "iDRW_prop(overlapped-GMP)")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_barplot(xlabs=xlabs, res_models, 0.6, perf_max+0.01)

############################### For performance comparison Line plot - Main figure ############################
# In iDRW_prop model, the model which has the best performance was selected in same Gamma. 

# load all models
gm_0.2 <- get(load('data/model/res_pa_GM_18_2_LOOCV.RData'))
gm_0.4 <- get(load('data/model/res_pa_GM_18_3_LOOCV.RData'))
gm_0.6 <- get(load('data/model/res_pa_GM_18_4_LOOCV.RData'))
gm_0.8 <- get(load('data/model/res_pa_GM_18_7_LOOCV.RData'))
gmr_0.2 <- get(load('data/model/res_pa_GMR_18_2_LOOCV.RData'))
gmr_0.4 <- get(load('data/model/res_pa_GMR_18_3_LOOCV.RData'))
gmr_0.6 <- get(load('data/model/res_pa_GMR_18_4_LOOCV.RData'))
gmr_0.8 <- get(load('data/model/res_pa_GMR_18_5_LOOCV.RData'))
gmr_d_0.2 <- get(load('data/model/res_pa_GMR_d_18_17_LOOCV.RData'))
gmr_d_0.4 <- get(load('data/model/res_pa_GMR_d_18_13_LOOCV.RData'))
gmr_d_0.6 <- get(load('data/model/res_pa_GMR_d_18_14_LOOCV.RData'))
gmr_d_0.8 <- get(load('data/model/res_pa_GMR_d_18_15_LOOCV.RData'))
gr_0.2 <- get(load('data/model/res_pa_GR_28_0.2_LOOCV.RData'))
gr_0.4 <- get(load('data/model/res_pa_GR_28_0.4_LOOCV.RData'))
gr_0.6 <- get(load('data/model/res_pa_GR_28_0.6_LOOCV.RData'))
gr_0.8 <- get(load('data/model/res_pa_GR_28_0.8_LOOCV.RData'))
gr_d_0.2 <- get(load('data/model/res_pa_GR_d_29_1_LOOCV.RData'))
gr_d_0.4 <- get(load('data/model/res_pa_GR_d_29_6_LOOCV.RData'))
gr_d_0.6 <- get(load('data/model/res_pa_GR_d_29_11_LOOCV.RData'))
gr_d_0.8 <- get(load('data/model/res_pa_GR_d_29_16_LOOCV.RData'))


# Plot GMR
title <- c("Performance comparison under varying Gamma")
xlabs <- rep(c("iDRW(GM)", "iDRW(GMP)", "iDRW_prop(GMP)",
               "iDRW(GP)", "iDRW_prop(GP)"), each=4)
Gamma_list <- rep(c("0.2", "0.4", "0.6", "0.8"), 5)

res_models <- list(gm_0.2, gm_0.4, gm_0.6, gm_0.8,
                   gmr_0.2, gmr_0.4, gmr_0.6, gmr_0.8,
                   gmr_d_0.2, gmr_d_0.4, gmr_d_0.6, gmr_d_0.8,
                   gr_0.2, gr_0.4, gr_0.6, gr_0.8,
                   gr_d_0.2, gr_d_0.4, gr_d_0.6, gr_d_0.8)

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_lineplot_multi(title, xlabs, res_models, perf_min-0.01, perf_max+0.01, Gamma_list)
#############################################################################################################


################## For performance comparison (result 18) - Main figure #######################################
# 730*472
gmr_mean <- get(load('data/model/res_pa_GMR_mean.RData'))
gmr_median <- get(load('data/model/res_pa_GMR_median.RData'))
gmr_concat <- get(load('data/model/res_pa_GMR_concat.RData'))
gm <- get(load('data/model/res_pa_GM_18_3_LOOCV.RData'))
gmr <- get(load('data/model/res_pa_GMR_18_3_LOOCV.RData'))
gmr_d <- get(load('data/model/res_pa_GMR_d_18_17_LOOCV.RData'))
gr <- get(load('data/model/res_pa_GR_28_0.6_LOOCV.RData'))
gr_d <- get(load('data/model/res_pa_GR_d_29_16_LOOCV.RData'))

res_models <- list(gmr_mean, gmr_median, gmr_concat, gm, gmr, gmr_d, gr, gr_d)

xlabs <- c("Mean", "Median", "Concat",
           "iDRW(GM)", "iDRW(GMP)", "iDRW_prop(GMP)",
           "iDRW(GP)", "iDRW_prop(GP)")

group <- c("{p(G),p(M),p(P)}", "{p(G),p(M),p(P)}", "{p(G),p(M),p(P)}",
           "{p(GM)}", "{p(GMP)}", "{p(GMP)}", "{p(GP)}", "{p(GP)}")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_barplot(xlabs = xlabs, res_models = res_models, perf_min = perf_min-0.01, perf_max = 0.75, group=group)

########## Performance of GP (result 28) ########
# 500*520
gmr <- get(load('data/model/res_pa_GMR_18_3_LOOCV.RData'))
gr <- get(load('data/model/res_pa_GR_28_0.6_LOOCV.RData'))
gr_d <- get(load('data/model/res_pa_GR_d_29_16_LOOCV.RData'))

res_models <- list(gmr, gr, gr_d)

xlabs <- c("iDRW(GMP)", "iDRW(GP)", "iDRW_prop(GP)")

perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
perf_barplot(xlabs, res_models, perf_min-0.01, 0.75)

#-------- confusion matrix, specificity, sensitivity
# GM
# Accuracy : 0.7021277
# sensitivity : 0.7085, specificity : 0.6949
# AUC : 0.701729
gm_pred <- gm$pred$pred[which(gm$pred$mtry == gm$bestTune$mtry)]
gm_obs <- gm$pred$obs[which(gm$pred$mtry == gm$bestTune$mtry)]
gm_pred <- as.numeric(as.character(gm_pred))
gm_obs <- as.numeric(as.character(gm_obs))

gm_conf <- confusionMatrix(table(gm_pred, gm_obs))

gm_AUC_perf <- performance(prediction(gm_pred, gm_obs), "auc")
gm_AUC <- gm_AUC_perf@y.values[[1]]

# GMR
# Accuracy : 0.7340426
# sensitivity : 0.7236, specificity : 0.7458
# AUC : 0.7346904
gmr_pred <- gmr$pred$pred[which(gmr$pred$mtry == gmr$bestTune$mtry)]
gmr_obs <- gmr$pred$obs[which(gmr$pred$mtry == gmr$bestTune$mtry)]
gmr_pred <- as.numeric(as.character(gmr_pred))
gmr_obs <- as.numeric(as.character(gmr_obs))

gmr_conf <- confusionMatrix(table(gmr_pred, gmr_obs))

gmr_AUC_perf <- performance(prediction(gmr_pred, gmr_obs), "auc")
gmr_AUC <- gmr_AUC_perf@y.values[[1]]


# GMR_d
# Accuracy : 0.7287234
# sensitivity : 0.7286, specificity : 0.7288
# AUC : 0.7287284
gmr_d_pred <- gmr_d$pred$pred[which(gmr_d$pred$mtry == gmr_d$bestTune$mtry)]
gmr_d_obs <- gmr_d$pred$obs[which(gmr_d$pred$mtry == gmr_d$bestTune$mtry)]
gmr_d_pred <- as.numeric(as.character(gmr_d_pred))
gmr_d_obs <- as.numeric(as.character(gmr_d_obs))

gmr_d_conf <- confusionMatrix(table(gmr_d_pred, gmr_d_obs))

gmr_d_AUC_perf <- performance(prediction(gmr_d_pred, gmr_d_obs), "auc")
gmr_d_AUC <- gmr_d_AUC_perf@y.values[[1]]
