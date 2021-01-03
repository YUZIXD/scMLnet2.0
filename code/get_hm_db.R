library(dplyr)
library(reshape2)
library(ggplot2)
library(cowplot)
options(stringsAsFactors = F)

######################
## prepare LigRecDB ##
######################

load("./prior_knowledge/LigRec/scMLnet.human.LigRec.rda")
LigRec.DB <- hm.LigRec
LigRec.DB$secondary.source[grep("ppi_prediction_go",LigRec.DB$secondary.source)] <- "ppi_prediction_go_nichenet"
LigRec.DB$secondary.source[grep("ppi_prediction$",LigRec.DB$secondary.source)] <- "ppi_prediction_nichenet"
LigRec.DB$secondary.source[grep("PMID|PMC",LigRec.DB$secondary.source)] <- "text_minning_cellchat"
LigRec.DB$secondary.source[grep("KEGG",LigRec.DB$secondary.source)] <- "KEGG_lr_cellchat"

LigRec.DB.sum <- as.data.frame.matrix(table(LigRec.DB$secondary.source, LigRec.DB$primary.source))
LigRec.DB.sum$data.source <- rownames(LigRec.DB.sum)
LigRec.DB.sum <- melt(LigRec.DB.sum,id.vars = "data.source")
LigRec.DB.sum <- LigRec.DB.sum[!LigRec.DB.sum$value==0,]
p1 <- ggplot(LigRec.DB.sum, aes(x = data.source,y = value,fill = variable))+
  geom_bar(stat ="identity")+ theme_bw() + 
  labs(x = "",y = "number of LR pairs", title = "Human LigRecDB")+  
  geom_text(aes(label = value),size = 5,vjust = -0.25)+ 
  guides(fill = guide_legend(reverse = F))+  
  theme(plot.title = element_text(size = 18, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(angle = -30,hjust = 0.1,vjust = 1,size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.position = 'bottom', 
        legend.key.size=unit(0.8,'cm')) 

LigRec.DB <- LigRec.DB %>% dplyr::select(source, target, support = secondary.source) 
rownames(LigRec.DB) <- NULL

#####################
## prepare RecTFDB ##
#####################

load("./rior_knowledge/RecTF/scMLnet.human.ppi.rda")
RecTF.DB <- hm.ppi
RecTF.DB <- RecTF.DB[!(RecTF.DB$src_type == "CHEBI" | RecTF.DB$dest_type == "CHEBI"),]
RecTF.DB <- RecTF.DB[!duplicated(RecTF.DB),]
table(RecTF.DB$type, RecTF.DB$direction)

# binding interaction / undirected interaction
add_pairs <- RecTF.DB[RecTF.DB$direction == "undirected" | RecTF.DB$type == "Process(binding/association)",] 
add_pairs <- add_pairs[,c(3:4,1:2,5:12)]
colnames(add_pairs) <- colnames(RecTF.DB)
RecTF.DB <- rbind(RecTF.DB, add_pairs)
RecTF.DB <- RecTF.DB[!duplicated(RecTF.DB),]

RecTF.DB <- RecTF.DB %>% dplyr::select(source = src, target = dest, primary.source = pathway, secondary.source = database)
RecTF.DB$secondary.source[grep("KEGG",RecTF.DB$secondary.source)] <- "KEGG_ppi_graphite"

RecTF.DB.sum <- as.data.frame(table(RecTF.DB$secondary.source))
p2 <- ggplot(RecTF.DB.sum, aes(x = Var1,y = Freq,fill = Var1))+
  geom_bar(stat ="identity")+ theme_bw() + 
  labs(x = "",y = "number of LR pairs", title = "Human RecTFDB")+ 
  geom_text(aes(label = Freq),size = 5,vjust = -0.25)+ 
  guides(fill = guide_legend(reverse = F))+ 
  theme(plot.title = element_text(size = 18, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(angle = -30,hjust = 0.1,vjust = 1,size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12),
        legend.position = 'bottom', 
        legend.key.size=unit(0.8,'cm')) 

RecTF.DB <- RecTF.DB %>% dplyr::select(source, target, support = secondary.source, pathway = primary.source)
rownames(RecTF.DB) <- NULL

####################
## prepare TFTGDB ##
####################

load("./prior_knowledge/TFTG/scMLnet.human.TFTG.rda")
TFTG.DB <- hm.TFTG
TFTG.DB$secondary.source[intersect(grep("PMID|PMC",TFTG.DB$secondary.source),grep("TRRUST",TFTG.DB$primary.source))] <- "text_minning_trrust"
TFTG.DB$secondary.source[intersect(grep("PMID|PMC",TFTG.DB$secondary.source),grep("HTRIdb",TFTG.DB$primary.source))] <- "text_minning_htridb"
TFTG.DB$secondary.source[grep("kegg",TFTG.DB$secondary.source)] <- "KEGG_tftg_Regnetwork"
TFTG.DB$secondary.source <- gsub("ensembl","Ensembl",TFTG.DB$secondary.source)
TFTG.DB$secondary.source <- gsub("hprd","HPRD",TFTG.DB$secondary.source)
TFTG.DB$secondary.source <- gsub("tred","TRED",TFTG.DB$secondary.source)
TFTG.DB$secondary.source <- gsub("ucsc","UCSC",TFTG.DB$secondary.source)
table(TFTG.DB$secondary.source, TFTG.DB$primary.source)

TFTG.DB.sum <- as.data.frame.matrix(table(TFTG.DB$secondary.source, TFTG.DB$primary.source))
TFTG.DB.sum$data.source <- rownames(TFTG.DB.sum)
TFTG.DB.sum <- melt(TFTG.DB.sum,id.vars = "data.source")
TFTG.DB.sum <- TFTG.DB.sum[!TFTG.DB.sum$value==0,]
p3 <- ggplot(TFTG.DB.sum, aes(x = data.source,y = value,fill = variable))+
  geom_bar(stat ="identity")+ theme_bw() + 
  labs(x = "",y = "number of TFTG pairs", title = "Human TFTGDB")+   
  geom_text(aes(label = value),size = 5,vjust = -0.25)+ 
  guides(fill = guide_legend(reverse = F))+  
  theme(plot.title = element_text(size = 18, vjust = 0.5, hjust = 0.5),
        axis.text.x = element_text(angle = -30,hjust = 0.1,vjust = 1,size = 12),
        axis.title = element_text(size = 12),
        legend.title = element_blank(), 
        legend.text = element_text(size = 12), 
        legend.position = 'bottom', 
        legend.key.size=unit(0.8,'cm')) 

TFTG.DB <- TFTG.DB %>% dplyr::select(source, target, support = secondary.source)
rownames(TFTG.DB) <- NULL

plot_grid(p1,p2,p3, ncol = 3)

###########
## check ##
###########

## check Receptors in ppi
Receptors <- LigRec.DB %>% dplyr::select(target) %>% unlist() %>% unique()
Receptors = Receptors[Receptors %in% unique(c(RecTF.DB$source, RecTF.DB$target))]

LigRec.DB = LigRec.DB %>% filter(target %in% Receptors)
Ligands <- LigRec.DB %>% dplyr::select(source) %>% unlist() %>% unique()

## check TFs in ppi
TFs <- TFTG.DB %>% dplyr::select(source) %>% unlist() %>% unique()
TFs = TFs[TFs %in% unique(c(RecTF.DB$source, RecTF.DB$target))]

TFTG.DB = TFTG.DB %>% filter(source %in% TFs)

save(Ligands,Receptors,TFs,file = "./prior_knowledge/prior_hm_nodes.rda")
save(LigRec.DB, RecTF.DB, TFTG.DB,file = "./prior_knowledge/prior_hm_edges.rda")

#############################
## obtain weighted network ##
#############################

## support function ##

#refer to 'construct_weighted_networks' function in NicheNet
creat_weighted_networks <- function(LigRec.DB, RecTF.DB, TFTG.DB, support_weights_df) {
  if (!is.data.frame(LigRec.DB)) 
    stop("LigRec.DB must be a data frame")
  if (!is.data.frame(RecTF.DB)) 
    stop("RecTF.DB must be a data frame")
  if (!is.data.frame(TFTG.DB)) 
    stop("TFTG.DB must be a data frame")
  if (!is.data.frame(support_weights_df) || sum((support_weights_df$weight > 1)) != 0) 
    stop("support_weights_df must be a data frame and no data support weight may be higher than 1")
   
  requireNamespace("dplyr")
  
  ## remove support for which weight equals 0
  support_weights_df = support_weights_df %>% filter(weight > 0)
  
  ## create weighted network
  LigRec_weight_network = LigRec.DB %>% 
    inner_join(support_weights_df, by = "support") %>% 
    group_by(source, target) %>% 
    summarize(weight = sum(weight)) %>% 
    ungroup()
  
  RecTF_weight_network = RecTF.DB %>% 
    inner_join(support_weights_df, by = "support") %>% 
    group_by(source, target) %>% 
    summarize(weight = sum(weight)) %>% 
    ungroup()
  
  TFTG_weight_network = TFTG.DB %>% 
    inner_join(support_weights_df, by = "support") %>% 
    group_by(source, target) %>% 
    summarize(weight = sum(weight)) %>% 
    ungroup()
  
  weighted_networks = list(LigRec = LigRec_weight_network, 
                           RecTF = RecTF_weight_network, 
                           TFTG = TFTG_weight_network)
  
  return(weighted_networks)
} 

#refer to 'apply_hub_corrections' function in NicheNet
mitigate_hub_effect <- function(weighted_networks, LigRec_hub, RecTF_hub, TFTG_hub){
  
  if (!is.list(weighted_networks)) 
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$LigRec)) 
    stop("LigRec must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$RecTF)) 
    stop("RecTF must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$TFTG)) 
    stop("TFTG must be a data frame or tibble object")
  if (!is.numeric(weighted_networks$LigRec$weight)) 
    stop("LigRec must contain a column named `weight`")
  if (!is.numeric(weighted_networks$RecTF$weight)) 
    stop("RecTF must contain a column named `weight`")
  if (!is.numeric(weighted_networks$TFTG$weight)) 
    stop("TFTG must contain a column named `weight`")
  
  if (LigRec_hub < 0 | LigRec_hub > 1) 
    stop("LigRec_hub must be a number between 0 and 1 (0 and 1 included)")
  if (RecTF_hub < 0 | RecTF_hub > 1) 
    stop("RecTF_hub must be a number between 0 and 1 (0 and 1 included)")
  if (TFTG_hub < 0 | TFTG_hub > 1) 
    stop("TFTG_hub must be a number between 0 and 1 (0 and 1 included)")
  
  requireNamespace("dplyr")
  
  ## get weight network
  LigRec_network = weighted_networks$LigRec
  RecTF_network = weighted_networks$RecTF
  TFTG_network = weighted_networks$TFTG
  
  ## mitigate the effect of hub node
  if (LigRec_hub > 0) {
    LigRec_network = LigRec_network %>% group_by(target) %>% count(target) %>% ungroup() %>% 
      inner_join(LigRec_network,., by = "target") %>% group_by(source) %>% 
      mutate(weight = weight/(n^LigRec_hub)) %>% 
      dplyr::select(-n) %>% 
      ungroup()
  }
  if (RecTF_hub > 0) {
    RecTF_network = RecTF_network %>% group_by(target) %>% count(target) %>% ungroup() %>% 
      inner_join(RecTF_network,., by = "target") %>% group_by(source) %>% 
      mutate(weight = weight/(n^RecTF_hub)) %>% 
      dplyr::select(-n) %>% 
      ungroup()
  }
  if (TFTG_hub > 0) {
    TFTG_network = TFTG_network %>% group_by(target) %>% count(target) %>% ungroup() %>% 
      inner_join(TFTG_network,., by = "target") %>% group_by(source) %>% 
      mutate(weight = weight/(n^TFTG_hub)) %>% 
      dplyr::select(-n) %>% 
      ungroup()
  }
  
  return(list(LigRec = LigRec_network, 
              RecTF = RecTF_network,
              TFTG = TFTG_network))
}

#get the convert network for calculation of potential matrix
convert_weighted_network <- function(weighted_networks, network){
  
  if (!is.list(weighted_networks)) 
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$LigRec)) 
    stop("LigRec must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$RecTF)) 
    stop("RecTF must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$TFTG)) 
    stop("TFTG must be a data frame or tibble object")
  if (!is.numeric(weighted_networks$LigRec$weight)) 
    stop("LigRec must contain a column named 'weight'")
  if (!is.numeric(weighted_networks$RecTF$weight)) 
    stop("RecTF must contain a column named 'weight'")
  if (!is.numeric(weighted_networks$TFTG$weight)) 
    stop("TFTG must contain a column named 'weight'")
  
  if (network != "LigRec" & network != "RecTF" & network != "TFTG") 
    stop("network must be 'LigRec' or 'RecTF' or 'TFTG'")
  
  requireNamespace("dplyr")
  
  ## get weight network
  LigRec_network = weighted_networks$LigRec
  RecTF_network = weighted_networks$RecTF
  TFTG_network = weighted_networks$TFTG
  
  ## convert ids to numeric
  allgenes = c(LigRec_network$source, LigRec_network$target,
               RecTF_network$source, RecTF_network$target,
               TFTG_network$source, TFTG_network$target) %>% unique() %>% sort()
  allgenes_integer = allgenes %>% factor() %>% as.numeric()
  allgenes_id_tbl = data.frame(allgenes, allgenes_integer) %>% as_tibble()
  
  mapper = function(df, value_col, name_col) setNames(df[[value_col]], df[[name_col]])
  id2allgenes = mapper(allgenes_id_tbl, "allgenes_integer", "allgenes")
  
  ## select weighted network
  weighted_network = weighted_networks[[network]]
  
  ## convert weighted network 
  weighted_network = weighted_network %>% 
    mutate(from_allgenes = id2allgenes[source], to_allgenes = id2allgenes[target]) %>% 
    arrange(from_allgenes) %>% 
    dplyr::select(from_allgenes, to_allgenes, weight)
  
  return(list(weighted_network = weighted_network,
              id2allgenes = id2allgenes))
}

#refer to 'construct_tf_target_matrix' function in NicheNet
construct_adjacency_matrix <- function (weighted_networks, network, source_as_cols = FALSE, 
                                        standalone_output = TRUE){
  
  if (!is.list(weighted_networks)) 
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$LigRec)) 
    stop("LigRec must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$RecTF)) 
    stop("RecTF must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$TFTG)) 
    stop("TFTG must be a data frame or tibble object")
  if (!is.numeric(weighted_networks$LigRec$weight)) 
    stop("LigRec must contain a column named 'weight'")
  if (!is.numeric(weighted_networks$RecTF$weight)) 
    stop("RecTF must contain a column named 'weight'")
  if (!is.numeric(weighted_networks$TFTG$weight)) 
    stop("TFTG must contain a column named 'weight'")
  
  if (network != "LigRec" & network != "RecTF" & network != "TFTG") 
    stop("network must be 'LigRec' or 'RecTF' or 'TFTG'")
  if (!is.logical(source_as_cols) | length(source_as_cols) != 1) 
    stop("source_as_cols must be a logical vector of length 1")
  if (!is.logical(standalone_output) | length(standalone_output) != 1) 
    stop("standalone_output must be a logical vector of length 1")
  
  requireNamespace("dplyr")
  
  ## get convert weighted network
  convert_result = convert_weighted_network(weighted_networks, network)
  weighted_network = convert_result$weighted_network
  allgenes = names(convert_result$id2allgenes)
  
  ## get adjacency matrix
  adjacency_matrix = Matrix::sparseMatrix(weighted_network$from_allgenes %>% as.integer, 
                                    weighted_network$to_allgenes %>% as.integer, 
                                    x = weighted_network$weight %>% as.numeric, 
                                    dims = c(length(allgenes), length(allgenes)))
  rownames(adjacency_matrix) = allgenes
  colnames(adjacency_matrix) = allgenes
  
  ## only keep nodes in specific weighted network
  if (standalone_output == TRUE) {
    source = weighted_networks[[network]] %>% .$source %>% unique()
    target = weighted_networks[[network]] %>% .$target %>% unique()
    adjacency_matrix = adjacency_matrix[source, target]
  }
  
  ## make sure source nodes as column
  if (source_as_cols == TRUE) {
    adjacency_matrix = adjacency_matrix %>% as.matrix()
    adjacency_matrix = adjacency_matrix %>% t()
  }
  
  return(adjacency_matrix)
} 

#refer to 'construct_rec_tf_matrix' function in NicheNet
RWR_wrapper = function(seed, G, delta, id2allgenes) {
  
  # always restart in specific seeds
  # prepare preference vector P(the initial probability distribution), only the seeds have values different from zero in P.
  P = rep(0, times = length(igraph::V(G)))
  P[id2allgenes[seed]] = 1
  
  # get a proximity measure from every graph node to the specific seed.
  partial_matrix = igraph::page_rank(G, algo = c("prpack"), vids = igraph::V(G),directed = TRUE, damping = delta, personalized = P) %>% .$vector
  rwr_matrix = matrix(unlist(partial_matrix), ncol = length(P), byrow = TRUE)
  
  return(rwr_matrix)
  
}
SPL_wrapper = function(seed, G, id2allgenes) {
  
  # calculate spl distance between seed and every other node in graph
  distances = igraph::distances(graph = G, v = id2allgenes[seed], to = igraph::V(G),mode = "out")
  
  # the shorter the better
  # change na and infinite predictions to -1 in order to later replace them by the maximum (= no probability)
  distances[is.nan(distances)] = -1
  distances[is.infinite(distances)] = -1
  
  distances[distances == -1] = max(distances)
  
  # reverse distances to weights: short distance: high weight
  spl_matrix = max(distances) - distances
  
  return(spl_matrix)
  
}
Direct_wrapper = function(seed, G, id2allgenes) {
  
  direct_matrix = G[id2allgenes[seed],]
  
  direct_matrix = matrix(direct_matrix, nrow = 1)
  
  return(direct_matrix)
}

#apply network algorithm to correct the RecTF potential matrix
construct_correct_matrix <- function(weighted_networks, seeds_of_interset, 
                                     algorithm = "RWR", restart_probability = 0.7, 
                                     source_as_cols = FALSE, standalone_output = TRUE) {
  
  if (!is.list(weighted_networks)) 
    stop("weighted_networks must be a list object")
  if (!is.data.frame(weighted_networks$LigRec)) 
    stop("LigRec must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$RecTF)) 
    stop("RecTF must be a data frame or tibble object")
  if (!is.data.frame(weighted_networks$TFTG)) 
    stop("TFTG must be a data frame or tibble object")
  if (!is.numeric(weighted_networks$LigRec$weight)) 
    stop("LigRec must contain a column named 'weight'")
  if (!is.numeric(weighted_networks$RecTF$weight)) 
    stop("RecTF must contain a column named 'weight'")
  if (!is.numeric(weighted_networks$TFTG$weight)) 
    stop("TFTG must contain a column named 'weight'")
  
  if (!is.vector(seeds_of_interset)) 
    stop("seeds_of_interset must be a vector object")
  if (sum((unique(seeds_of_interset) %in% unique(c(weighted_networks$RecTF$source, weighted_networks$RecTF$target))) == FALSE) > 0) 
    warning("One or more seeds of interest not present in the RecTF weight network.")
  
  if (algorithm != "RWR" & algorithm != "SPL" & algorithm != "direct") 
    stop("algorithm must be 'RWR' or 'SPL' or 'direct'")
  if (algorithm == "RWR") {
    if (restart_probability < 0 | restart_probability >= 1) 
      stop("restart_probability must be a number between 0 and 1 (0 included, 1 not)")
  }
  
  if (!is.logical(source_as_cols) | length(source_as_cols) != 1) 
    stop("source_as_cols must be a logical vector of length 1")
  if (!is.logical(standalone_output) | length(standalone_output) != 1) 
    stop("standalone_output must be a logical vector of length 1")
  
  requireNamespace("dplyr")
  
  ## get convert weighted network
  convert_result = convert_weighted_network(weighted_networks, network = "RecTF")
  RecTF_network = convert_result$weighted_network
  id2allgenes = convert_result$id2allgenes
  allgenes = names(id2allgenes)

  ## apply network algorithm
  if (algorithm == "RWR") {
    
    ## construct adjacency matrix 
    RecTF_matrix = Matrix::sparseMatrix(RecTF_network$from_allgenes %>% as.integer, 
                                        RecTF_network$to_allgenes %>% as.integer, 
                                        x = RecTF_network$weight %>% as.numeric, 
                                        dims = c(length(allgenes), length(allgenes)))
    
    ## construct weight graph from adjacency matrix
    RecTF_igraph = igraph::graph_from_adjacency_matrix(RecTF_matrix, 
                                                       weighted = TRUE, 
                                                       mode = "directed")
    
    # restart in every seeds individual
    complete_matrix = lapply(seeds_of_interset, RWR_wrapper, RecTF_igraph, 
                             restart_probability, id2allgenes)
    
  }
  else if (algorithm == "SPL") {
    
    # the greater the weight, the smaller the distance, the greater the potential
    RecTF_network = RecTF_network %>% mutate(weight = 1/weight)
    
    # construct adjacency matrix 
    RecTF_matrix = Matrix::sparseMatrix(RecTF_network$from_allgenes %>% as.integer, 
                                        RecTF_network$to_allgenes %>% as.integer, 
                                        x = RecTF_network$weight %>% as.numeric, 
                                        dims = c(length(allgenes), length(allgenes)))
    
    # construct weight graph from adjacency matrix
    RecTF_igraph = igraph::graph_from_adjacency_matrix(RecTF_matrix, 
                                                       weighted = TRUE, 
                                                       mode = "directed")
    
    # get the proximity between the seed and all the other nodes in the graph.
    complete_matrix = lapply(seeds_of_interset, SPL_wrapper, RecTF_igraph, id2allgenes)
    
  }
  else if (algorithm == "direct") {
    
    # construct the network that link seed to seed
    seed_seed_network = tibble::tibble(from_allgenes = id2allgenes[seeds_of_interset], 
                                       to_allgenes = id2allgenes[seeds_of_interset]) %>% inner_join(RecTF_network %>% 
                                                                                                    filter(from_allgenes %in% id2allgenes[seeds_of_interset]) %>% 
                                                                                                    group_by(from_allgenes) %>% 
                                                                                                    top_n(1, weight) %>% 
                                                                                                    ungroup() %>% 
                                                                                                    distinct(from_allgenes, weight), by = "from_allgenes")
    # combine the RecTF_network and the network that link seed to seed
    RecTF_network = RecTF_network %>% bind_rows(seed_seed_network)
    
    # construct adjacency matrix 
    RecTF_matrix = Matrix::sparseMatrix(RecTF_network$from_allgenes %>% as.integer, 
                                        RecTF_network$to_allgenes %>% as.integer, 
                                        x = RecTF_network$weight %>% as.numeric, 
                                        dims = c(length(allgenes), length(allgenes)))
    
    # construct weight graph from adjacency matrix
    RecTF_igraph = igraph::graph_from_adjacency_matrix(RecTF_matrix, 
                                                       weighted = TRUE, 
                                                       mode = "directed")
    
    # get adjacency matrix of seeds_of_interset to other nodes in graph
    complete_matrix = lapply(seeds_of_interset, Direct_wrapper, RecTF_matrix, id2allgenes)
    
  }
  
  # merge 
  ReCTF_correct_matrix = do.call(rbind, complete_matrix)
  
  rownames(ReCTF_correct_matrix) = seeds_of_interset
  colnames(ReCTF_correct_matrix) = allgenes
  
  ## make sure source nodes as column
  if (source_as_cols == TRUE) {
    ReCTF_correct_matrix = ReCTF_correct_matrix %>% as.matrix()
    ReCTF_correct_matrix = ReCTF_correct_matrix %>% t()
  }
  
  ## only keep nodes in specific weighted network
  if (standalone_output == TRUE) {
    target = weighted_networks$TFTG %>% .$source %>% unique()
    ReCTF_correct_matrix = ReCTF_correct_matrix[, target]
  }
  
  return(as(ReCTF_correct_matrix, "sparseMatrix"))
}

## main ##

## load interaction
rm(list = ls())
load(file = "./prior_knowledge/prior_hm_edges.rda")

## support weight
support_weights_df <- data.frame(support=unique(c(LigRec.DB$support,
                                                  RecTF.DB$support,
                                                  TFTG.DB$support)), 
                                 weight=1)

## parameter used for optimization
hyperparameter_list = list(LigRec_hub = 1, # mitigation of the effect of hub node in LigRec weighted network
                           RecTF_hub = 1, # mitigation of the effect of hub node in RecTF weighted network
                           TFTG_hub = 1, # mitigation of the effect of hub node in TFTG weighted network
                           restart_probability = 0.7) # parameter in Random walk with restart algorithms

## create weighted network
weighted_networks = creat_weighted_networks(LigRec.DB, RecTF.DB, TFTG.DB, support_weights_df)

## mitigate the importance of hubs node
weighted_networks = mitigate_hub_effect(weighted_networks, 
                                        LigRec_hub = hyperparameter_list$LigRec_hub, 
                                        RecTF_hub = hyperparameter_list$RecTF_hub, 
                                        TFTG_hub = hyperparameter_list$TFTG_hub)

## construct adjacency matrix of LigRec and TFTG weight network
LigRec_matrix = construct_adjacency_matrix(weighted_networks, network = "LigRec")
TFTG_matrix = construct_adjacency_matrix(weighted_networks, network = "TFTG")

## construct adjacency matrix of RecTF weight network
Receptors <- weighted_networks$LigRec %>% select(target) %>% unlist() %>% unique()
RecTF_matrix_RWR = construct_correct_matrix(weighted_networks, seeds_of_interset = Receptors, 
                                            algorithm = "RWR")

## get potential matrices list
potential_matrixs = list(LigRec_matrix = LigRec_matrix,
                         RecTF_matrix = RecTF_matrix_RWR,
                         TFTG_matrix = TFTG_matrix)

saveRDS(potential_matrixs,file = "./prior_knowledge/hm_potential_matrix_RWR.rds")
saveRDS(weighted_networks,file = "./prior_knowledge/hm_weighted_networks.rds")
