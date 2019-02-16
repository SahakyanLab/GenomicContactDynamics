################################################################################
#Density plot to show distribution of the number of bases per bin with age data

################################################################################
#functions to source
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
#objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
#HiC_Human21 persistent contact files
#persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
persist.dir = "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
#BaseAgeVsPersist files
#baseAge.persist.dir = paste0(objective.dir, "/out_BaseAgeVsPersist_test")
baseAge.persist.dir = paste0(objective.dir, "/out_BaseAgeVsPersist")
dir.create(paste0(objective.dir, "/out_AgeCountPerNtis"))
output.dir = paste0(objective.dir, "/out_AgeCountPerNtis")  

library(foreach)
#library(ggpubr) #also loads ggplot2
library(ggplot2)

#2(2MB gap) or "05"(0.5 MB gap), refers to minimum gap accepted to classify 
#a contact, two points should be far enough to filter for contacts within a TAD
gcb.v = c("2", "05")

chr.v = c(1:22, "X")

plotOnly <- FALSE 

################################################################################
foreach(gcb=gcb.v, .inorder=TRUE) %do% {
  for(chr in chr.v){
      load(paste0(persist.dir, "/chr", chr, "_Persist_min", gcb, "Mb.RData"))
      mx <- cbind(i=as.numeric(PERSIST.MX$hits[,"i"]),
                  j=as.numeric(PERSIST.MX$hits[,"j"]),
                  ntis=as.numeric(PERSIST.MX$ntis))
      rm( list=c("PERSIST.MX") )
      gc()
      
      #take unique bins per ntis
      #bins.uniq.ntis <- by( data=mx[, c("i", "j")], INDICES=mx[,"ntis"], 
      #                     FUN=function(x){unique(unlist(x))} )
      #ntis.uniq <- as.numeric(names(bins.uniq.ntis))
      
      #take bins per ntis
      bins.ntis <- by( data=mx[, c("i", "j")], INDICES=mx[,"ntis"], 
                       FUN=unlist )
      ntis.uniq <- as.numeric(names(bins.ntis))
      
      load(paste0(baseAge.persist.dir, "/min", gcb, "Mb_chr", chr,  
                  "_bin40Kb_BaseAgeVsPersist.RData"))
      #matrix for plotting
      lst <- lapply(X=as.list(ntis.uniq), FUN=function(x){
        perNtis <- bins.ntis[[x]]
        mtched <- match(perNtis, BASEAGE.BIN.MX[,c("bin")])
        ind.notNa <- which(!is.na(mtched))
        count.ntis <- BASEAGE.BIN.MX[mtched[ind.notNa],"count"]
        cbind(ntis=rep(as.numeric(x), length(count.ntis)), 
              bin=as.numeric(perNtis[ind.notNa]), 
              count=as.numeric(count.ntis))
      })
      rm( list=c("BASEAGE.BIN.MX") )
      gc()
      AGECOUNT.NTIS.MX <- do.call(rbind, lst)
      save(AGECOUNT.NTIS.MX, file=paste0(output.dir,"/plot_min", gcb, "Mb_chr", 
                                         chr, "_AgeCountPerNtis.RData"))
  }
}

#rm(list=ls())

