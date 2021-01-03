library(dplyr)
library(tidyr)
library(Hmisc)
library(org.Hs.eg.db)

######################
## prepare TRRUSTDB ##
######################

hm.trrust <- read.table("./prior_knowledge/TFTG/human_TRRUST.tsv",header = F,sep = "\t")
hm.trrust$primary.source <- rep("TRRUST",nrow(hm.trrust))
hm.trrust$V4 <- paste0("PMID:",hm.trrust$V4)
hm.trrust <- hm.trrust[,c(1:2,5,4)]
colnames(hm.trrust) <- c("source","target","primary.source","secondary.source")

######################
## prepare HTRIdbDB ##
######################

hm.htridb <- read.table("./prior_knowledge/TFTG/human_HTRIdb_TFTG.txt",
                        header = F,skip = 1,sep = "\t")
hm.htridb <- hm.htridb[,seq(1,14,by = 2)]
colnames(hm.htridb) <- read.table("./prior_knowledge/TFTG/human_HTRIdb_TFTG.txt",
                                  header = F,nrows = 1,sep = "\t")[-8]
hm.htridb$primary.source <- rep("HTRIdb",nrow(hm.htridb))
hm.htridb$PUBMED_ID <- paste0("PMID:",hm.htridb$PUBMED_ID)
hm.htridb <- hm.htridb[,c(3,5,8,7)]
colnames(hm.htridb) <- c("source","target","primary.source","secondary.source")

##########################
## prepare RegNetworkDB ##
##########################

hm.regnet <- read.csv("./prior_knowledge/TFTG/human_RegNetwork.csv")
table(hm.regnet$database)
hm.regnet <- hm.regnet[-grep("mir",hm.regnet$database),]
hm.regnet$primary.source <- rep("RegNetwork",nrow(hm.regnet))
hm.regnet <- hm.regnet[,c(1,3,8,5)]
colnames(hm.regnet) <- c("source","target","primary.source","secondary.source")
hm.regnet <- hm.regnet %>% mutate(secondary.source = strsplit(as.character(secondary.source), ",")) %>% unnest(secondary.source)
table(hm.regnet$secondary.source)

################
################

hm.TFTG <- rbind(rbind(hm.trrust,hm.htridb),hm.regnet)
save(hm.TFTG, file = "./prior_knowledge/TFTG/scMLnet.human.TFTG.rda")
