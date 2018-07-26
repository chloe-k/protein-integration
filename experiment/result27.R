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
# 

id_list <- c("27_0.2_G", "27_0.2_M", "27_0.2_R", "27_0.2_GM", "27_0.2_GP", "27_0.2_MP", "27_0.2_GMP")
type_list <- c("g", "m", "p", "gm", "gp", "mp", "gmp")

pack <- c("KEGGgraph", "igraph", "ggplot2", "annotate", "annotate", "org.Hs.eg.db", "diffusr", "DESeq2", "Matrix",
          "stringr", "caret", "e1071", "randomForest", "KEGG.db", "KEGGREST")


res_gmr_27_0.2 <- foreach(i=1:length(id_list), .packages = pack) %dopar%{
  make_GMR_model(id=id_list[i], prob = 0.001, Gamma = 0.2, type_used = type_list[i])
}