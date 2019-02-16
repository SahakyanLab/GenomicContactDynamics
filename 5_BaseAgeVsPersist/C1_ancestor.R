#BaseAgeAncestorVsPersist
start_time <- Sys.time()

###############################################################################  
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for task
#objective.dir <- "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
#HiC_Human21 persistent contact files
#persist.dir = "/Users/ltamon/Database/HiC_features_GSE87112_RAWpc"
persist.dir = "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
#BaseAgeVsPersist files
#baseAge.persist.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist/out_BaseAgeVsPersist_test"
baseAge.persist.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist/out_BaseAgeVsPersist"
#dir.create(paste0(objective.dir, "/out_ancestor_test"))
#output.dir = paste0(objective.dir, "/out_ancestor_test") 
#dir.create(paste0(objective.dir, "/out_ancestor_noOutliers"))
output.dir = paste0(objective.dir, "/out_ancestor") 

library(ggplot2)
library(foreach)
#install.packages(doParallel)
#install.packages(iterators)
library(itertools) #isplitVector

nCPU=4

#2(2Mb) or "05"(0-5 Mb), minimum gap acceptedbetween contacting bins 
gc.v = "2" #c("2", "05") 

chr.v = c("21")

#minimum allowed number of age data per bin
minAgeCountPerBin <- 1L 

plotOnly <- FALSE

gc="2"
chr="21"

################################################################################   
source(paste0(lib, "/UTL_NCPU.R"))
source(paste0(lib, "/findHPM.R"))
source(paste0(lib, "/ggLayers_persist.R"))
#source(paste0(lib, "/multiplot.R"))

#give data points per box in boxplot
#cite
give.n <- function(x){ return(c(y=-0.1, label = length(x))) }

myplot <- function(df=df){
  ggplot(data=as.data.frame(df), aes(x=ntis, y=df[,grp])) +
    #geom_boxplot(outlier.size = NA, aes(fill=factor(ntis)), na.rm=TRUE) +
    geom_boxplot(outlier.size = 0.1, aes(fill=factor(ntis)), na.rm=TRUE) +
    guides(fill=FALSE) + 
    #stat_summary(geom="text", label=totalperntis,
    #             fun.y=function(x){max(x, na.rm)}, vjust=-2, size=2) + 
    stat_summary(fun.data=give.n, geom="text", size=2) + 
    ggtitle(paste0("chr", chr, " _min", gc, "Mb_", grp)) +
    scale_x_discrete(name="Number (non-0 Tissues/Cell lines)", 
                     limits=ntis.uniq) +
    scale_fill_manual(values=PersistScoreColour(ntis.uniq)) +
    scale_y_continuous(name="Fraction of Human/Primate/Mammal-specific base",
                       limits=c(-0.1, 1.1)) +
    theme.persist #from sourced "ggLayers_persist.R"
}

################################################################################  
################################################################################  
for(gc in gc.v){
  for(chr in chr.v){
    
    if(plotOnly==FALSE){
      
      load(file=paste0(baseAge.persist.dir, "/min", gc, "Mb_chr", chr, 
                       "_bin40Kb_BaseAgeVsPersist.RData"))
      #filter for bins that satisfy minAgeCountPerBin
      incl.bin <- BASEAGE.BIN.MX[BASEAGE.BIN.MX[,"count"]>=minAgeCountPerBin,"bin"]
      load(paste0(persist.dir, "/chr", chr, "_Persist_min", gc, "Mb.RData"))
      #initialize matrix for fraction of HPM, same order with PERSIST.MX$hits
      HPM.NTIS.MX <- cbind(i=as.numeric(PERSIST.MX$hits[,"i"]),
                           j=as.numeric(PERSIST.MX$hits[,"j"]),
                           HumanSpecificBase=rep(NA),
                           PrimateSpecificBase=rep(NA),
                           MammalSpecificBase=rep(NA),
                           TotalNumberOfBase=rep(NA),
                           ntis=PERSIST.MX$ntis)
      ntis.uniq <- sort( as.numeric(unique(PERSIST.MX$ntis)) )
      rm( list=c("PERSIST.MX") )
      gc()
      incl.ij.ind <- as.numeric( which(HPM.NTIS.MX[,"i"]%in%incl.bin & HPM.NTIS.MX[,"j"]%in%incl.bin) )
      
      #add fraction of HPM and total base per contact to HPM.NTIS.MX
      incl.ij.ind.len <- length(incl.ij.ind)
      
      #### FOREACH EXECUTION #########
      findHPM.out <- foreach( itr=isplitVector(1:incl.ij.ind.len, chunks=nCPU), 
               .combine="rbind", .inorder=TRUE,
               .export=c("BASEAGE.BIN.MX", 
                         "HPM.NTIS.MX", 
                         "incl.ij.ind"),
               .noexport=ls()[!ls()%in%c("BASEAGE.BIN.MX", 
                                         "HPM.NTIS.MX", 
                                         "incl.ij.ind")]) 
      %op% {
        #combine ancestor of i and j as a string and calculate fraction
        findHPM.chunk <- sapply(X=itr, simplify=FALSE, FUN=function(x){
          ancestor.str <- paste(BASEAGE.BIN.MX[ as.character(HPM.NTIS.MX[incl.ij.ind[x],"i"]), 
                                                "ancestor" ],
                                BASEAGE.BIN.MX[ as.character(HPM.NTIS.MX[incl.ij.ind[x],"j"]), 
                                                "ancestor" ], 
                                sep=";")
          findHPM(ancestor.str)
        })
        
        return(do.call(rbind, findHPM.chunk))
      }
      ### END OF FOREACH EXECUTION ###
      
      #add the calculated fractions to the initialized matrix
      HPM.NTIS.MX[incl.ij.ind, 3:6] <- findHPM.out
      
      save(HPM.NTIS.MX, 
           file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_", 
                       minAgeCountPerBin, "_HPMfraction_AncestorVsPersist.RData"))
      
    } else {
      load(file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_", 
                       minAgeCountPerBin, "_HPMfraction_AncestorVsPersist.RData"))
    }
    
    #calculate total number of base age data per ntis
    #totalperntis <- by( data=HPM.NTIS.MX[,"TotalNumberOfBase"], 
    #                    INDICES=HPM.NTIS.MX[,"ntis"],
    #                    FUN=sum, simplify=TRUE )
    #ntis.uniq <- as.numeric(names(totalperntis))
    #totalperntis <- as.vector(totalperntis)
    
    #-------------------------------------------
    #myplots.list <- list()
    for(grp in c("HumanSpecificBase", "PrimateSpecificBase", "MammalSpecificBase")){
      df <- HPM.NTIS.MX[,c(grp, "ntis")]
      p <- myplot(df=df)
      #myplots.list[[as.character(grp)]] <- p
      ggsave(file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_",
                         minAgeCountPerBin, "_", grp, 
                         "_AncestorVsPersist.jpeg"),
             plot=p)
      
    }
    #pdf( file=paste0(output.dir, "/plot_min", gc, "Mb_chr", chr, "_",
    #                  minAgeCountPerBin, "_HPM_AncestorVsPersist.pdf"),
    #     width=13, height=8.5 )
    #multiplot(plotlist=myplots.list, cols=3)
    #dev.off()
  }
}
#rm(list=ls())
end_time <- Sys.time()
end_time-start_time   

