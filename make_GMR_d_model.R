make_GMR_d_model <- function(id){
  # id - 18_28
  
  msg <- paste(c("id is : ",id), collapse = '')
  print(msg)
  
  sapply(file.path("utils",list.files("utils", pattern="*.R")),source)
  
  # make directory
  datapath <- file.path('data')
  if(!dir.exists(datapath)) dir.create(datapath)
  
  gdacpath <- file.path(datapath, 'BRCA_GDAC')
  
  respath <- file.path('result')
  if(!dir.exists(respath)) dir.create(respath)
  
  
  # read RNAseq, Methylation data, RPPA data
  #data_all_path <- file.path(datapath, "data.RData")
  data_all_path <- file.path(datapath, "Entrez_data.RData")
  if(!file.exists(data_all_path)) {
    year <- 3
    read_data(year, datapath)
  }
  load(data_all_path)
  
  
  # read rda data
  # graphpath <- file.path(datapath,'directGraph.rda')
  # pathSetpath <- file.path(datapath,'pathSet.rda')
  graphpath <- file.path(datapath,'directGraph(Entrez).rda')
  pathSetpath <- file.path(datapath,'pathSet(Entrez).rda')
  if(!(file.exists(graphpath) && file.exists(pathSetpath))) {
    print('directGraph and pathSet do not exist, now creating directGraph and pathSet is start')
    cons_KEGGgraph(datapath)
  }
  load(file.path(graphpath))
  load(file.path(pathSetpath))
  
  
  
  dppipath <- file.path(datapath, 'DppiGraph(Entrez).rda')  # directed edge PPI
  # dppipath <- file.path(datapath, 'DppiGraph_rdc.rda')  # directed edge PPI
  # dppipath <- file.path(datapath, 'DppiGraph_W_str.rda')
  
  if(!file.exists(dppipath)) {
    print('ppiGraph does not exist, now creating ppiGraph start')
    cons_ppi(datapath, gdacpath, rppa)
  }
  load(file.path(dppipath))
  
  
  # directed pathway graph provided in DRWPClass
  g <- directGraph 
  V(g)$name <- paste("g",V(g)$name,sep="")
  
  m <- directGraph
  V(m)$name <-paste("m",V(m)$name,sep="")
  
  r <- directGraph
  V(r)$name <-paste("p",V(r)$name,sep="")
  
  p <- DppiGraph
  V(p)$name <-paste("p",V(p)$name,sep="")
  
  y=list(good_samples, poor_samples)
  
  #------------------------- RNAseq + Methyl + RPPA(PPI Graph) -------------------------#
  gmr <- list(g, m, r)
  testStatistic <- c("DESeq2", "t-test", "t-test")
  profile_name <- c("rna(Entrez)", "meth(Entrez)", "rppa(diffused_Pathway_Graph_Entrez)")
  x=list(rnaseq, imputed_methyl, rppa)
  
  # model_name -> res_pa_GMR_d_18_28.RData
  # id -> result18_28_GMR_d
  result_name <- paste(c('result',id,'_GMR_d'), collapse = '')
  
  # fit.iDRWPClass(x=x, y=y, globalGraph=gmp, testStatistic= testStatistic, profile_name = profile_name,
  #                datapath = datapath, respath = respath, pathSet=pathSet, method = "DRW", samples = samples,
  #                id = "result18_28_GMR_d", prob = 0.8, Gamma = 0.4, pranking = "t-test", mode = "GMP", AntiCorr=FALSE, DEBUG=TRUE)
  
  model <- fit.classification(y=y, samples = samples, id = result_name, datapath = datapath, respath = respath,
                              profile_name = profile_name, method = "DRW", pranking = "t-test", classifier = "rf",
                              nFolds = 5, numTops=50, iter = 10)
  
  
  model_path <- paste(c('data/model/res_pa_GMR_d_',id,'.RData'), collapse = '')
  
  name <- paste(c('res_pa_GMR_d_', id), collapse='')
  assign(x = name, value = model)
  
  save(list=name, file=file.path(model_path))
  msg <- paste(c(result_name,' is done'), collapse = '')
  print(msg)
  
  return(model)
}
