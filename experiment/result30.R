# parameter tuning for result25_5

################################## Result 30 ############################################################
num_cores <- detectCores()/2
registerDoParallel(cores = num_cores)


id_list <- c("30_0.2", "30_0.4", "30_0.6", "30_0.8", "30_0.9")
lim_list <- c(5)
Gamma_list <- c(0.2, 0.4, 0.6, 0.8, 0.9)


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")
res_gmr_25 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], lim=lim_list[1], prob = 0.001, Gamma = Gamma_list)
}

res_models <- list()
for(i in 1:length(id_list)){
  model <- get(load(paste(c('data/model/res_pa_GMR_', id_list[i], '.RData'), collapse = '')))
  res_models <- c(res_models, list(model))
}

# Plot GMR
title <- c("Result 30 GMR")
xlabs <- c("0.2", "0.4", "0.6", "0.8", "0.9")

perf_min <- min(sapply(X = res_models, FUN = function(x){mean(x$resample$Accuracy)}))
perf_max <- max(sapply(X = res_models, FUN = function(x){mean(x$resample$Accuracy)}))
perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.01, perf_max = perf_max+0.01)
# perf_lineplot(fname_res = 'result/res_loocv_tuneK.txt', perf_min=55, perf_max=95)