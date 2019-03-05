# Set minimum age count per contacting bin

###############################################################################
#main directory for task
#objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
#BaseAgeVsPersist files
baseAge.persist.dir =  paste0(objective.dir, "/out_BaseAgeVsPersist")
#HPM fractions per contact files
hpm.ntis.dir = paste0(objective.dir, "/out_ancestor") 

library(foreach)
library(doParallel) 

nCPU=23
registerDoParallel(cores=nCPU)

#2(2Mb) or "05"(0.5 Mb), minimum gap acceptedbetween contacting bins
gcb = "05"

chr.v = c(1:22, "X")
#chr = "1"

minAgeCountPerBin = 3000L #3000L

###############################################################################
###############################################################################
foreach(chr=chr.v, .inorder = FALSE,
        .export=c("gcb", 
                  "minAgeCountPerBin",
                  "objective.dir",
                  "baseAge.persist.dir",
                  "hpm.ntis.dir"),
        .noexport=ls()[!ls()%in%c("gcb", 
                                  "minAgeCountPerBin",
                                  "objective.dir",
                                  "baseAge.persist.dir",
                                  "hpm.ntis.dir")] 
        ) %dopar% {
          
  load(file=paste0(baseAge.persist.dir, "/min", gcb, "Mb_chr", chr, 
                   "_bin40Kb_BaseAgeVsPersist.RData"))
  #filter for bins that satisfy minAgeCountPerBin
  incl.bin <- BASEAGE.BIN.MX[BASEAGE.BIN.MX[,"count"]>=minAgeCountPerBin,"bin"]
  rm( list=c("BASEAGE.BIN.MX") )
  gc()
  load(file=paste0(hpm.ntis.dir, "/plot_min", gcb, "Mb_chr", chr, 
                   "_1_HPMfraction_AncestorVsPersist.RData"))
  #contacts with HPM fraction information (minAgeCountPerBin=1)
  notNA.ij.ind <- which(!is.na(HPM.NTIS.MX[,"TotalNumberOfBase"]))
  incl.ij.ind <- which(HPM.NTIS.MX[,"i"]%in%incl.bin & HPM.NTIS.MX[,"j"]%in%incl.bin)
  changetoNA.ind <- as.numeric( setdiff(notNA.ij.ind, incl.ij.ind) )
  HPM.NTIS.MX[changetoNA.ind, 3:6] <- rep(NA, 4)
  
  save(HPM.NTIS.MX, 
       file=paste0(hpm.ntis.dir, "/plot_min", gcb, "Mb_chr", chr, "_", 
                   minAgeCountPerBin, "_HPMfraction_AncestorVsPersist.RData"))
}

    

  


 
