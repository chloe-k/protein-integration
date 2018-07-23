################################## Result 25 ############################################################

num_cores <- detectCores()/3
num_cores <- floor(num_cores*2)
registerDoParallel(cores = num_cores)


id_list <- c("25_5", "25_10", "25_15")
lim_list <- c(5, 10, 15)


pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")
res_gmr_25 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], lim=lim_list[i])
}

res_gmr_d_25 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_d_model(id=id_list[i], lim = lim_list[i])
}




#########################################################################################################################################
# Plot for Result 25

# res_models <- list(res_pa_GMR_25_5_LOOCV, res_pa_GMR_d_25_5_LOOCV, 
#                    res_pa_GMR_25_10_LOOCV, res_pa_GMR_d_25_10_LOOCV,
#                    res_pa_GMR_25_15_LOOCV, res_pa_GMR_d_25_15_LOOCV)
# 
# title <- c("Result 25")
# xlabs <- c("GMR_5", "GMR_d_5", "GMR_10", "GMR_d_10", "GMR_15", "GMR_d_15")
# perf_min <- min(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
# perf_max <- max(sapply(X = res_models, FUN = function(x){max(x$results$Accuracy)}))
# perf_boxplot(title, xlabs, res_models, perf_min = perf_min-0.15, perf_max = perf_max+0.15)

res_model_25_5 <- list(res_pa_GMR_25_5_LOOCV, res_pa_GMR_d_25_5_LOOCV)
res_model_25_10 <- list(res_pa_GMR_25_10_LOOCV, res_pa_GMR_d_25_10_LOOCV)
res_model_25_15 <- list(res_pa_GMR_25_15_LOOCV, res_pa_GMR_d_25_15_LOOCV)

title <- c("Result 25_5")
xlabs <- c("GMR", "GMR_d")
perf_min <- min(sapply(X = res_model_25_5, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_5, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_5, perf_min = perf_min-0.15, perf_max = perf_max+0.15)


title <- c("Result 25_10")
xlabs <- c("GMR", "GMR_d")
perf_min <- min(sapply(X = res_model_25_10, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_10, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_10, perf_min = perf_min-0.15, perf_max = perf_max+0.15)


title <- c("Result 25_15")
xlabs <- c("GMR", "GMR_d")
perf_min <- min(sapply(X = res_model_25_15, FUN = function(x){max(x$results$Accuracy)}))
perf_max <- max(sapply(X = res_model_25_15, FUN = function(x){max(x$results$Accuracy)}))
perf_boxplot(title, xlabs, res_model_25_15, perf_min = perf_min-0.15, perf_max = perf_max+0.15)
###############################################################
# GMR model

# make_GMR_model <- function(id, lim){
#   # id - 18_28
#   
#   msg <- paste(c("id is : ",id), collapse = '')
#   print(msg)
#   
#   sapply(file.path("utils",list.files("utils", pattern="*.R")),source)
#   
#   # make directory
#   datapath <- file.path('data')
#   if(!dir.exists(datapath)) dir.create(datapath)
#   
#   gdacpath <- file.path(datapath, 'BRCA_GDAC')
#   
#   respath <- file.path('result')
#   if(!dir.exists(respath)) dir.create(respath)
#   
#   
#   # read RNAseq, Methylation data, RPPA data
#   #data_all_path <- file.path(datapath, "data.RData")
#   data_all_path <- file.path(datapath, "Entrez_data.RData")
#   if(!file.exists(data_all_path)) {
#     year <- 3
#     read_data(year, datapath)
#   }
#   load(data_all_path)
#   
#   
#   # read rda data
#   # graphpath <- file.path(datapath,'directGraph.rda')
#   # pathSetpath <- file.path(datapath,'pathSet.rda')
#   graphpath <- file.path(datapath,'directGraph(Entrez).rda')
#   pathSetpath <- file.path(datapath,'pathSet(Entrez).rda')
#   if(!(file.exists(graphpath) && file.exists(pathSetpath))) {
#     print('directGraph and pathSet do not exist, now creating directGraph and pathSet is start')
#     cons_KEGGgraph(datapath)
#   }
#   load(file.path(graphpath))
#   load(file.path(pathSetpath))
#   
#   
#   
#   dppipath <- file.path(datapath, 'DppiGraph(Entrez).rda')  # directed edge PPI
#   # dppipath <- file.path(datapath, 'DppiGraph_rdc.rda')  # directed edge PPI
#   # dppipath <- file.path(datapath, 'DppiGraph_W_str.rda')
#   
#   if(!file.exists(dppipath)) {
#     print('ppiGraph does not exist, now creating ppiGraph start')
#     cons_ppi(datapath, gdacpath, rppa)
#   }
#   load(file.path(dppipath))
#   
#   
#   # directed pathway graph provided in DRWPClass
#   g <- directGraph 
#   V(g)$name <- paste("g",V(g)$name,sep="")
#   
#   m <- directGraph
#   V(m)$name <-paste("m",V(m)$name,sep="")
#   
#   r <- directGraph
#   V(r)$name <-paste("p",V(r)$name,sep="")
#   
#   p <- DppiGraph
#   V(p)$name <-paste("p",V(p)$name,sep="")
#   
#   y=list(good_samples, poor_samples)
#   
#   #------------------------- RNAseq + Methyl + RPPA(Pathway Graph) -------------------------#
#   gmr <- g %du% m %du% r
#   testStatistic <- c("DESeq2", "t-test", "t-test")
#   profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(Pathway_Graph_Entrez)")
#   x=list(rnaseq, imputed_methyl, rppa)
#   
#   # model_name -> res_pa_GMR_18_28.RData
#   # id -> result18_28_GMR
#   result_name <- paste(c('result',id,'_GMR'), collapse = '')
#   
#   fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                  datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, lim = lim,
#                  id = result_name, prob = 0.001, Gamma = 0.4, pranking = "t-test", mode = "GMR", AntiCorr=FALSE, DEBUG=TRUE)
#   
#   model <- fit.classification(y=y, samples = samples, id = result_name, datapath = datapath, respath = respath,
#                               profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                               nFolds = 5, numTops=50, iter = 10)
#   
#   
#   model_path <- paste(c('data/model/res_pa_GMR_',id,'_LOOCV.RData'), collapse = '')
#   
#   name <- paste(c('res_pa_GMR_', id, '_LOOCV'), collapse='')
#   assign(x = name, value = model)
#   
#   save(list=name, file=file.path(model_path))
#   msg <- paste(c(result_name,' is done'), collapse = '')
#   print(msg)
#   
#   return(model)
# }



###############################################################
# GMR_d model

# make_GMR_d_model <- function(id, lim){
#   # id - 18_28
#   
#   msg <- paste(c("id is : ",id), collapse = '')
#   print(msg)
#   
#   sapply(file.path("utils",list.files("utils", pattern="*.R")),source)
#   
#   # make directory
#   datapath <- file.path('data')
#   if(!dir.exists(datapath)) dir.create(datapath)
#   
#   gdacpath <- file.path(datapath, 'BRCA_GDAC')
#   
#   respath <- file.path('result')
#   if(!dir.exists(respath)) dir.create(respath)
#   
#   
#   # read RNAseq, Methylation data, RPPA data
#   #data_all_path <- file.path(datapath, "data.RData")
#   data_all_path <- file.path(datapath, "Entrez_data.RData")
#   if(!file.exists(data_all_path)) {
#     year <- 3
#     read_data(year, datapath)
#   }
#   load(data_all_path)
#   
#   
#   # read rda data
#   # graphpath <- file.path(datapath,'directGraph.rda')
#   # pathSetpath <- file.path(datapath,'pathSet.rda')
#   graphpath <- file.path(datapath,'directGraph(Entrez).rda')
#   pathSetpath <- file.path(datapath,'pathSet(Entrez).rda')
#   if(!(file.exists(graphpath) && file.exists(pathSetpath))) {
#     print('directGraph and pathSet do not exist, now creating directGraph and pathSet is start')
#     cons_KEGGgraph(datapath)
#   }
#   load(file.path(graphpath))
#   load(file.path(pathSetpath))
#   
#   
#   
#   dppipath <- file.path(datapath, 'DppiGraph(Entrez).rda')  # directed edge PPI
#   # dppipath <- file.path(datapath, 'DppiGraph_rdc.rda')  # directed edge PPI
#   # dppipath <- file.path(datapath, 'DppiGraph_W_str.rda')
#   
#   if(!file.exists(dppipath)) {
#     print('ppiGraph does not exist, now creating ppiGraph start')
#     cons_ppi(datapath, gdacpath, rppa)
#   }
#   load(file.path(dppipath))
#   
#   
#   # directed pathway graph provided in DRWPClass
#   g <- directGraph 
#   V(g)$name <- paste("g",V(g)$name,sep="")
#   
#   m <- directGraph
#   V(m)$name <-paste("m",V(m)$name,sep="")
#   
#   r <- directGraph
#   V(r)$name <-paste("p",V(r)$name,sep="")
#   
#   p <- DppiGraph
#   V(p)$name <-paste("p",V(p)$name,sep="")
#   
#   y=list(good_samples, poor_samples)
#   
#   #------------------------- RNAseq + Methyl + RPPA(diffused Pathway Graph) -------------------------#
#   gmr <- list(g, m, r)
#   testStatistic <- c("DESeq2", "t-test", "t-test")
#   profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
#   x=list(rnaseq, imputed_methyl, rppa)
#   
#   # model_name -> res_pa_GMR_d_18_28.RData
#   # id -> result18_28_GMR_d
#   result_name <- paste(c('result',id,'_GMR_d'), collapse = '')
#   
#   fit.iDRWPClass(x=x, y=y, globalGraph=gmr, testStatistic= testStatistic, profile_name = profile_name,
#                  datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples, lim = lim,
#                  id = result_name, prob = 0.4, Gamma = 0.2, pranking = "t-test", mode = "GMR_d", AntiCorr=FALSE, DEBUG=TRUE)
#   
#   model <- fit.classification(y=y, samples = samples, id = result_name, datapath = datapath, respath = respath,
#                               profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
#                               nFolds = 5, numTops=50, iter = 10)
#   
#   
#   model_path <- paste(c('data/model/res_pa_GMR_d_',id,'_LOOCV.RData'), collapse = '')
#   
#   name <- paste(c('res_pa_GMR_d_', id, '_LOOCV'), collapse='')
#   assign(x = name, value = model)
#   
#   save(list=name, file=file.path(model_path))
#   msg <- paste(c(result_name,' is done'), collapse = '')
#   print(msg)
#   
#   return(model)
# }


