cons_graph <- function(datapath){
  
  # This code implements constructing global directed gene-gene graph from all kegg pathway(525 pathways)
  # 525 pathways -> 327 pathways(It is the number of existing kgml)
  # directGraph
  # node : 7389
  # edge : 58399 
  
  datapath <- file.path('data')
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
    cat("kgml files already exists")
    kgml_list <- list.files(path=kegg_db, full.names=TRUE)
  }
  
  
  graphs <- list()
  pathSet <- list()
  
  for(kgml in kgml_list){
    
    # parse kgml file
    k <- parseKGML(kgml)
    
    # pathway set
    node_names <- unlist(lapply(X=nodes(k),
                                FUN=function(x) {
                                  ifelse(getType(x)=="gene", getSYMBOL(sapply( strsplit(getName(x),":"), '[[', 2 ), data='org.Hs.eg'), NA)
                                }), use.names = F)
    pathSet[[substring(getName(k),9)]] <- node_names[!is.na(node_names)]
    
    # graph set
    pathwayG <- parseKGML2Graph(kgml, expandGenes=TRUE)
    graphs <- append(graphs, pathwayG)
  }
  
  # make direct graph
  merged_G <- mergeGraphs(graphs)
  directGraph <- igraph.from.graphNEL(merged_G)
  
  V(directGraph)$name <- sapply(X = V(directGraph)$name, FUN = function(x) getSYMBOL(strsplit(x, ":")[[1]][2], data='org.Hs.eg'), USE.NAMES = F)
  
  save(directGraph, file=file.path(datapath, paste(c("directGraph", "rda"), collapse='.')))
  save(pathSet, file=file.path(datapath, paste(c("pathSet", "rda"), collapse='.')))

}

