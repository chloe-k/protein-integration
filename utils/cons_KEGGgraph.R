cons_KEGGgraph <- function(datapath){
  
  # This code implements constructing global directed gene-gene graph from all kegg pathway(525 pathways)
  # 525 pathways -> 327 pathways(It is the number of existing kgml)
  # directGraph
  # node : 7389
  # edge : 58399 
  
  kegg_db <- file.path(datapath,'KEGG_DB')
  
  if(!dir.exists(kegg_db)) {
    cat("kgml file doesn't exist.")
    kegg_list <- scan("../data/kegg_pathway", what="", sep = "\n")
    kgml_list <- list()
    
    for(pathid in kegg_list){
      tryCatch(
        {
          tmp <- paste('../data/KEGG_DB/hsa',pathid,sep="")
          tmp <- paste(tmp,".xml", sep="")
          retrieveKGML(pathwayid = pathid, organism = "hsa", destfile = tmp, method = "wget")
          kgml_list <- append(kgml_list, tmp)
          cat(pathid,'is downloaded!...................')
          
        },
        error=function(e){
          cat(pathid,'is not exist')
          message(e)
        }
        
      )
    }
    
  } else{
    cat("kgml files already exists\n")
    kgml_list <- list.files(path='data/KEGG_DB', full.names=TRUE)
  }
  
  graphs <- list()
  pathSet<-list()
  
  for(kgml in kgml_list){
    # parse kgml file
    k <- parseKGML(kgml)
    
    # pathway set
    node_names <- c()
    for(node in nodes(k)) {
      if(getType(node) == "gene") {
        #node_names <- c(node_names, getSYMBOL(sapply( strsplit(getName(node),":"), '[[', 2 ), data='org.Hs.eg'))
        node_names <- sapply( strsplit(getName(node),":"), '[[', 2 )
      }
    }
    
    pathSet[[substring(getName(k),9)]] <- unique(node_names)
    
    # graph set
    pathwayG <- KEGGpathway2Graph(k, expandGenes=TRUE)
    graphs <- append(graphs, pathwayG)
  }
  
  merged_G <- mergeGraphs(graphs)
  directGraph <- igraph.from.graphNEL(merged_G)
  
  #V(directGraph)$name <- sapply(X = V(directGraph)$name, FUN = function(x) getSYMBOL(strsplit(x, ":")[[1]][2], data='org.Hs.eg'), USE.NAMES = F)
  V(directGraph)$name <- sapply(X = V(directGraph)$name, FUN = function(x) strsplit(x, ":")[[1]][2], USE.NAMES = F)
  
  # save(directGraph, file=file.path(datapath, paste(c("directGraph", "rda"), collapse='.')))
  # save(pathSet, file=file.path(datapath, paste(c("pathSet", "rda"), collapse='.')))
  
  save(directGraph, file=file.path(datapath, paste(c("directGraph(Entrez)", "rda"), collapse='.')))
  save(pathSet, file=file.path(datapath, paste(c("pathSet(Entrez)", "rda"), collapse='.')))
  
}

