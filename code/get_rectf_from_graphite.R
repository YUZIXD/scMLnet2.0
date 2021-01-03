if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("graphite")
BiocManager::install("org.Hs.eg.db")

library(graphite) #October 27, 2020
library(dplyr)
library(org.Hs.eg.db)

#########################################
## get pathwayDB from graphite package ##
#########################################

## check databases and species
pathwayDatabases()
human.dbs <- pathwayDatabases() %>% dplyr::filter(species == "hsapiens") %>% dplyr::select(database) %>% unlist()

## run in parallel
options(Ncpus = 3)

## get human interaction
t3 <- Sys.time()
hm.dbs <- lapply(human.dbs, function(h.db){
  pws <- pathways(species = "hsapiens", database = h.db)
  pws.sym <- convertIdentifiers(pws, "SYMBOL")
  
  pws.dat <- lapply(pws.sym, function(pw.sym){
    # pw.sym <- pws.sym[[2]]
    pw.dat <- edges(pw.sym)
    if(nrow(pw.dat)){
      pw.dat <- pw.dat %>% mutate(numnodes = pw.dat %>% dplyr::select(src,dest) %>% unlist() %>% unique() %>% length(),
                                  numedges = nrow(pw.dat),
                                  id = pw.sym@id, 
                                  pathway = pw.sym@title, 
                                  database = pw.sym@database,
                                  species = pw.sym@species)
    }
    pw.dat
  })
  pws.dat <- do.call(rbind,pws.dat)
  pws.dat
})
hm.ppi <- do.call(rbind,hm.dbs)
t4 <- Sys.time()
t4-t3
save(hm.ppi, file = "./prior_knowledge/RecTF/scMLnet.human.ppi.rda")
