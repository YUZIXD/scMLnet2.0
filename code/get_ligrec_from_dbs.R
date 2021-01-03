library(dplyr)
library(tidyr)
library(Hmisc)

########################
## prepare cellchatDB ##
########################

load("./prior_knowledge/LigRec/human_CellChatDB.rda")
hm.cellchat <- CellChatDB.human$interaction
hm.cellchat <- hm.cellchat[,c("ligand","receptor","evidence")]

hm.cellchat.complex <- CellChatDB.human$complex
hm.cellchat.complex$combine <- apply(hm.cellchat.complex,1,function(x){
  x = x[nchar(x)!=0]
  paste(x,collapse = "_")
})

hm.cellchat$ligand <- lapply(hm.cellchat$ligand,function(x){
  x = ifelse(x %in% rownames(hm.cellchat.complex),hm.cellchat.complex[x,5],x)
})
hm.cellchat$receptor <- lapply(hm.cellchat$receptor,function(x){
  x = ifelse(x %in% rownames(hm.cellchat.complex),hm.cellchat.complex[x,5],x)
})

hm.cellchat <- hm.cellchat %>% mutate(ligand = strsplit(as.character(ligand), "_")) %>% unnest(ligand)
hm.cellchat <- hm.cellchat %>% mutate(receptor = strsplit(as.character(receptor), "_")) %>% unnest(receptor)

hm.cellchat$primary.source <- rep("CellChat",nrow(hm.cellchat))
hm.cellchat <- hm.cellchat %>% dplyr::select(source = ligand, target = receptor, primary.source,secondary.source = evidence)

##########################
## prepare connectomeDB ##
##########################

hm.coondb <- read.csv("E:/study_r/prior_knowledge/LigRec/human_mouse_connectomeDB2020.csv",header = T)
hm.coondb$primary.source <- rep("connectomeDB2020",nrow(hm.coondb))
hm.coondb <- hm.coondb[,c(4,7,11,2)]
colnames(hm.coondb) <- c("source","target","primary.source","secondary.source")
hm.coondb$secondary.source <- unlist(lapply(hm.coondb$secondary.source, function(x){
  switch(x,
         `Ramilowski_2015_Literature_supported` = "Ramilowski_2015",
         `Hou et al. 2020 (this publication)` = "Hou_2020",           
         `Baccin et al. 2020 (RNA-Magnet)` = "RNA-Magnet",             
         `Efremova et al. 2020 (CellphoneDB)` = "CellphoneDB",       
         `Cabello-Aguilar et al. 2020 (SingleCellSignalR)` = "SingleCellSignalR",
         `No?l et al. 2020 (ICELLNET)` = "ICELLNET")
}))
table(hm.coondb$secondary.source)

########################
## prepare NicheNetDB ##
########################

hm.nichenet <- read.csv("E:/study_r/prior_knowledge/LigRec/human_mouse_NicheNet.csv")
table(hm.nichenet$database)
hm.nichenet <- hm.nichenet[hm.nichenet$database %in% c("ppi_prediction","ppi_prediction_go"),]
hm.nichenet$primary.source <- rep("NicheNet",nrow(hm.nichenet))
hm.nichenet <- hm.nichenet[,c(1,2,5,4)]
colnames(hm.nichenet) <- c("source","target","primary.source","secondary.source")

###############
###############

hm.LigRec <- rbind(rbind(hm.nichenet,hm.coondb),hm.cellchat)
save(hm.LigRec, file = "./prior_knowledge/LigRec/scMLnet.human.LigRec.rda")



