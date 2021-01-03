library(Matrix)
library(dplyr)
library(Seurat)

####################
## supportting function ##
####################

# get barcode list of cluster of interset and other cluster
getBarList <- function(Aclu, GCMat, BarCluTable)
{
  AcluBar <- BarCluTable %>% filter(.,Cluster == Aclu) %>% select(.,Barcode) %>% unlist()
  names(AcluBar) <- NULL
  
  AllBar <- BarCluTable %>% select(Barcode) %>% unlist()
  OtherBar <- setdiff(AllBar,AcluBar)
  
  result <- list(AcluBar,OtherBar)
  return(result)
}

# perform fisher test to get activate TFs or pathways
fisher_test <- function(subset1,subset2,backgrond)
{
  a=length(intersect(subset1,subset2))
  b=length(subset1)-a
  c=length(subset2)-a
  d=length(backgrond)-a-b-c
  matrix=matrix(c(a,c,b,d),nrow=2)
  fisher.test(matrix,alternative="greater")$p.value
}

# get nodes in prior Database
getNodeList <- function(Database, Nodetype)
{
  
  NodeList <- Database[,Nodetype] %>% unlist() %>% unique()
  
  return(NodeList)
  
}

##################
##    main function    ##
##################

runNormalize <- function(GCMat, norm.method = "LogNormalize")
{
  
  if (!(norm.method %in% c("LogNormalize", "CLR", "RC", "SCTransform"))) 
    stop("wrong normalization method!")
  
  seur = CreateSeuratObject(counts = GCMat)
  
  if(norm.method == "SCTransform") {
    seur = SCTransform(seur, verbose = TRUE)
    GCMat = seur@assays$SCT@data
  }else {
    seur = NormalizeData(seur, normalization.method = norm.method, verbose = TRUE)
    GCMat = seur@assays$RNA@data
  }
  
  return(GCMat)
}

getDiffExpGene <- function(GCMat, BarCluTable, Clus, min_pct = 0.1, min_diff.pct = -Inf,
                           logfc.cutoff = 0.25, p_val.cutoff = 0.05, Test.methods = "t")
{
  
  if (!is.logical(Raw.data) | length(Raw.data) != 1) 
    stop("Raw.data must be a logical vector of length 1")
  if (p_val.cutoff < 0 | p_val.cutoff > 1) 
    stop("p_val.cutoff must be a number between 0 and 1 (0 and 1 included)")
  if (!(Test.methods %in% c("wilcox", "bimod", "roc", "t", "MAST", "LR"))) 
    stop("wrong test method!")
  
  ## parameters
  cat(paste0("min_pct = ",min_pct,"\n"))
  cat(paste0("logfc.cutoff = ",logfc.cutoff,"\n"))
  cat(paste0("p_val.cutoff = ",p_val.cutoff,"\n"))
  cat(paste0("methods = ",Test.methods,"\n"))
  
  ## get barcode
  BarListResult <- getBarList(Aclu = Clus, GCMat, BarCluTable)
  Clus.1 <- BarListResult[[1]]
  Clus.2 <- BarListResult[[2]]
  
  ## find DEGs(use all other cells for FindMarkers)
  cat(paste0("get differentially expressed genes in ",Clus,"\n"))
  DEGs <- FindMarkers(object = GCMat, cells.1 = Clus.1, cells.2 = Clus.2, logfc.threshold = logfc.cutoff,
                      min.pct = min_pct, min.diff.pct = min_diff.pct, test.use = Test.methods, verbose = T)
  
  ## output
  up_gene <- DEGs %>% filter(p_val_adj <= p_val.cutoff & avg_logFC > 0) %>% rownames()
  #cat(paste0("find up-regulated genes: ",length(up_gene),"\n"))
  down_gene <- DEGs %>% filter(p_val_adj <= p_val.cutoff & avg_logFC < 0) %>% rownames()
  #cat(paste0("find down-regulated genes: ",length(down_gene),"\n"))
  
  return(DEGs)
  
}

getLigRec <- function(LigRec.DB, source_up, target_up)
{
  
  if (!is.data.frame(LigRec.DB)) 
    stop("LigRec.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(LigRec.DB)) 
    stop("LigRec.DB must contain a column named 'source'")
  if (!"target" %in% colnames(LigRec.DB)) 
    stop("LigRec.DB must contain a column named 'target'")
  
  # get ligand and receptor list
  LigGene <- LigRec.DB %>% select(source) %>% unlist() %>% unique()
  RecGene <- LigRec.DB %>% select(target) %>% unlist() %>% unique()
  TotLigRec <- paste(LigRec.DB$source, LigRec.DB$target, sep = "_") %>% unique()
  
  # get high expressed ligand and receptor
  LigHighGene <- intersect(LigGene,source_up)
  RecHighGene <- intersect(RecGene,target_up)
  
  # get activated LR pairs
  LRList <- paste(rep(LigHighGene,each = length(RecHighGene)),RecHighGene,sep = "_")
  LRList <- intersect(LRList,TotLigRec)
  
  # check result
  if(length(LRList)==0)
    stop("Error: No significant LigRec pairs")
  
  # get result
  LRTable <- LRList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(LRTable) <- c("source","target")
  
  cat(paste0("get ",length(LRList)," activated LR pairs\n"))
  return(LRTable)
  
}

getTFTG <- function(TFTG.DB, target.degs, target.genes)
{
  
  if (!is.data.frame(TFTG.DB)) 
    stop("TFTG.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(TFTG.DB)) 
    stop("TFTG.DB must contain a column named 'source'")
  if (!"target" %in% colnames(TFTG.DB)) 
    stop("TFTG.DB must contain a column named 'target'")
  
  # get TF list
  TF.list <- TFTG.DB %>% select(source) %>% unlist() %>% unique()
  
  # get Target list
  TG.list <- lapply(TF.list, function(x){
    TFTG.DB %>% filter(source == x)  %>% select(target) %>% unlist() %>% unique()
  })
  names(TG.list) <- TF.list
  
  # get target differently expressed genes
  DEGs <- target.degs
  
  # perform fisher test
  TFs <- lapply(TG.list, function(x){fisher_test(subset1 = x, subset2 = DEGs, backgrond = target.genes)})
  TFs <- unlist(TFs)
  TFs <- names(TFs)[TFs <= 0.05]
  TFs <- TFs[TFs %in% target.genes]
  
  # get activated LR pairs
  TFTGList <- TG.list[TFs]
  TFTGList <- lapply(TFTGList, function(x){intersect(x, DEGs)})
  TFTGList <- paste(rep(TFs, times = lengths(TFTGList)), unlist(TFTGList), sep = "_")
  
  # check result
  if(length(TFTGList)==0)
    stop("Error: No significant TFTG pairs")
  
  # get result
  TFTGTable <- TFTGList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(TFTGTable) <- c("source","target")
  
  cat(paste0("get ",length(TFTGList)," activated TFTG pairs\n"))
  return(TFTGTable)
  
}

getRecTF <- function(RecTF.DB, Rec.list, TF.list)
{
  
  if (!is.data.frame(RecTF.DB)) 
    stop("RecTF.DB must be a data frame or tibble object")
  if (!"source" %in% colnames(RecTF.DB)) 
    stop("RecTF.DB must contain a column named 'source'")
  if (!"target" %in% colnames(RecTF.DB)) 
    stop("RecTF.DB must contain a column named 'target'")
  
  # make sure Rec.list in RecTF.DB
  Rec.list <- Rec.list[Rec.list %in% RecTF.DB$source]
  Rec.list <- as.vector(Rec.list)
  
  # get TF activated by Receptors
  TFofRec <- lapply(Rec.list, function(x){
    RecTF.DB %>% filter(source == x)  %>% select(target) %>% unlist() %>% unique()
  })
  names(TFofRec) <- Rec.list
  
  # get all TF
  TFofALL <- RecTF.DB %>% select(target) %>% unlist() %>% unique()
  
  # perform fisher test
  Recs <- lapply(TFofRec, function(x){fisher_test(subset1 = x, subset2 = TF.list, backgrond = TFofALL)})
  Recs <- unlist(Recs)
  Recs <- names(Recs)[Recs <= 0.05]
  Recs <- Recs[Recs %in% target_gene]
  
  # get activated RecTF pairs
  RecTFList <- TFofRec[Recs]
  RecTFList <- lapply(RecTFList, function(x){intersect(x, TF.list)})
  RecTFList <- paste(rep(Recs, times = lengths(RecTFList)), unlist(RecTFList), sep = "_")
  
  # check result
  if(length(RecTFList)==0)
    stop("Error: No significant RecTF pairs")
  
  # get result
  RecTFTable <- RecTFList %>% strsplit(.,split = "_") %>% do.call(rbind, .) %>% as.data.frame()
  colnames(RecTFTable) <- c("source","target")
  
  cat(paste0("get ",length(RecTFList)," activated RecTF pairs\n"))
  return(RecTFTable)
  
}

################
##         input         ##
################

## input ##
data <- readRDS("./example/data.rds")
dim(data)
BarCluTable <- read.table("./example/barcodetype.txt", header = T, sep = "\t", stringsAsFactors = F)
table(BarCluTable$Cluster)

## parameters ##
LigClu <- "Mast cells"
RecClu <- "Secretory"

Raw.data = TRUE
min.cells = 3
min.features = 200
Norm.method <- "LogNormalize" #"CLR", "RC", "SCTransform"

# getDiffExpGene para
min_pct = 0.05
logfc.cutoff = 0.15 # gene select
p_val.cutoff = 0.05
Test.methods <- "t" #"bimod", "roc", "t", "MAST", "LR", "wilcox"

# get pairs para
all_mean <- rowMeans(GCMat)
hist(log10(all_mean), breaks=100, main="", col="grey80",
     xlab=expression(Log[10]~"average count"))
abline(v=log10(0.01), col="red", lwd=2, lty=2)
abundant.cutoff = 0.01
abs_logFC = 1

###################
##   prior knowledge   ##
###################

# user can input their own interaction DB which required two column: source and target
Databases <- readRDS("./prior_knowledge/hm_potential_matrix_RWR.rds")

LigRec.DB <- Databases$LigRec_matrix %>% as.matrix() %>% reshape2::melt() %>% 
  .[.$value > 0,] %>%
  select(source=Var1, target=Var2, score=value) %>% 
  mutate(source=as.vector(source), target=as.vector(target)) %>% as_tibble()

RecTF.DB <- Databases$RecTF_matrix %>% as.matrix() %>% reshape2::melt() %>% 
  .[.$value > quantile(.$value,0.9),] %>%
  select(source=Var1, target=Var2, score=value) %>% 
  mutate(source=as.vector(source), target=as.vector(target)) %>% as_tibble()

TFTG.DB <- Databases$TFTG_matrix %>% as.matrix() %>% reshape2::melt() %>%
  .[.$value > 0,] %>%
  select(source=Var1, target=Var2, score=value) %>% 
  mutate(source=as.vector(source), target=as.vector(target)) %>% as_tibble()

###############
##         run         ##
###############

## check ##
if (class(data)[1] == "matrix" & class(data)[1] == "dgCMatrix") 
  stop("data must be a matrix or sparse matrix")
if (!is.data.frame(BarCluTable)) 
  stop("BarCluTable must be a data frame or tibble object")
if (!is.character(BarCluTable$Barcode)) 
  stop("BarCluTable must contain a column named 'Barcode'")
if (!is.character(BarCluTable$Cluster)) 
  stop("BarCluTable must contain a column named 'Cluster'")
if (!(LigClu %in% unique(BarCluTable$Cluster)) & !(RecClu %in% unique(BarCluTable$Cluster)))
  stop("BarCluTable and RecClu must be a character in Cluster colunm of BarCluTable")

## keep the same barcodes order 
allCell <- intersect(colnames(data), BarCluTable$Barcode)
data <- data[,allCell]
BarCluTable <- BarCluTable[match(allCell,BarCluTable$Barcode),]

## filter ##
# filter genes on the number of cells expressing (eg.3)
if (min.cells > 0) {
  num.cells <- Matrix::rowSums(x = data > 0)
  data <- data[which(x = num.cells >= min.cells), ]
}
# Filter cells based on min.features (eg.200)
if (min.features > 0) {
  nfeatures <- Matrix::colSums(x = data > 0)
  data <- data[, which(x = nfeatures >= min.features)]
}

# update BarCluTable
BarCluTable <- BarCluTable[match(colnames(data),BarCluTable$Barcode),]

## normalization ##
if(Raw.data){
  cat("perform normalization\n")
  GCMat = runNormalize(data, norm.method = Norm.method)
}else{
  GCMat = data
}

## get DEGs of source cells and target cells ##
source_deg_tab <- getDiffExpGene(GCMat, BarCluTable, Clus = LigClu, min_pct = min_pct,
                                 logfc.cutoff = logfc.cutoff, Test.methods = Test.methods)
target_deg_tab <- getDiffExpGene(GCMat, BarCluTable, Clus = RecClu, min_pct = min_pct,
                                 logfc.cutoff = logfc.cutoff, Test.methods = Test.methods)

## get LigRec pairs ##
source_up <- source_deg_tab %>% filter(p_val_adj <= 0.05 & avg_logFC > 0) %>% rownames()
target_abundant <- names(all_mean)[all_mean > abundant.cutoff]
LigRecTab <- getLigRec(LigRec.DB, source_up, target_abundant)

## get TFTG pairs ##
target_gene <- rownames(GCMat)
target_deg <- target_deg_tab %>% filter(p_val_adj <= 0.05 & abs(avg_logFC) > abs_logFC) %>% rownames()
TFTGTab <- getTFTG(TFTG.DB, target_deg, target_gene)

## get RecTF pairs ##
Rec.list <- getNodeList(LigRecTab, "target")
TF.list <- getNodeList(TFTGTab, "source")
RecTFTab <- getRecTF(RecTF.DB, Rec.list, TF.list)

## updata ##
Receptors_in_Tab <- getNodeList(RecTFTab, "source")
LigRecTab <- LigRecTab[LigRecTab[,2] %in% Receptors_in_Tab,]

TFs_in_Tab <- getNodeList(RecTFTab, "target")
TFTGTab <- TFTGTab[TFTGTab[,1] %in% TFs_in_Tab,]

## output
result <- list("LigRec" = LigRecTab,
               "RecTF" = RecTFTab,
               "TFTar" = TFTGTab)

cat("Final LR pairs:",nrow(LigRecTab))
cat("Final RecTF pairs:",nrow(RecTFTab))
cat("Final TFTG pairs:",nrow(TFTGTab))

#########################
##        run end      ##
#########################