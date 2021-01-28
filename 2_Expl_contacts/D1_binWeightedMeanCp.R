################################################################################
# Calculate weighted mean Cp per bin, with and without considering Cp=0 (these 
# are long-range contacts not present in the cell)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/2_Expl_contacts"
    unmapbinPath = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation/out_unmappable/min2Mb_whole_Cs.norm_unmappable.csv"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/2_Expl_contacts"
    unmapbinPath = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation/out_unmappable/min2Mb_whole_Cs.norm_unmappable.csv"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
out.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binWeightedMeanCp")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
bin.len = 4e4
nCPU = 5L # Number of bins
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(gcb=="min2Mb"){
  mingap <- 50
  print(paste0("mingap=", mingap, " bins"), quote=FALSE)
} else if(gcb=="min05Mb"){
  mingap <- 12.5
  print(paste0("mingap=", mingap, " bins"), quote=FALSE)
} else {
  stop("Invalid input for gcb.")
}

unmapbin.df <- read.csv(file=unmapbinPath, header=TRUE, stringsAsFactors=FALSE)
chrLen.df <- read.delim(file=chrLenfile, header=TRUE, stringsAsFactors=FALSE)

#BINWMEANCP.DF <- list()
for(chr in chr.v){
  
  unmapbin <- unmapbin.df$unmapbin[unmapbin.df$chr==chr]
  unmapbin <- as.numeric(strsplit(x=unmapbin, split=";")[[1]])
  
  totbin <- ceiling(chrLen.df$length.bp[chrLen.df$chromosome==chr]/bin.len)
  if( totbin!=chrLen.df$bins.40kb[chrLen.df$chromosome==chr] ){
    stop("Checkpoint 1.")
  }
  
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  ij.df <- cbind.data.frame(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, Cp=PERSIST.MX$ntis)
  rm(PERSIST.MX); gc()
  
  toExport <- c("chr", "ij.df", "totbin", "mingap", "unmapbin")
  
  #### PARALLEL EXECUTION #########
  full <- foreach(itr=isplitVector(1:totbin, chunks=nCPU),
                  .inorder=TRUE, .combine="rbind",
                  .export=toExport, .noexport=ls()[!ls()%in%toExport]
                           
  ) %op% {
    
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(bin){
      
      v <- rep(x=0, times=totbin)
      names(v) <- 1:totbin
      
      # FILTER OUT bins in centromeric areas
      if(bin%in%unmapbin){
        
        v[v==0] <- NA
        return(v)
        
      } else {
        
        partn.i.TF <- ij.df$j==bin 
        partn.j.TF <- ij.df$i==bin 
        
        partn.i <- as.numeric(ij.df$i[partn.i.TF])
        partn.j <- as.numeric(ij.df$j[partn.j.TF])
        
        if( length(intersect(partn.i, partn.j))!=0 ){
          stop("Checkpoint 2.")
        }
        
        v[partn.i] <- ij.df$Cp[partn.i.TF]
        v[partn.j] <- ij.df$Cp[partn.j.TF]
        v[bin] <- NA
        
        rm(partn.j, partn.i, partn.j.TF, partn.i.TF)
        
        # FILTER OUT short-range contacts
        SR.bin <- 1:totbin%in%((bin-mingap):(bin+mingap))
        v[SR.bin] <- NA; rm(SR.bin)
        
        return(v)
        
      }
      
    })
    
    return( do.call("rbind", chunk) )
    
  }
  ### END OF PARALLEL EXECUTION ###
  
  rm(ij.df, unmapbin)
  
  wmeanCp0 <- rowMeans(x=full, na.rm=TRUE)
  full[full==0 & !is.na(full)] <- NA
  wmeanCp <- rowMeans(x=full, na.rm=TRUE)
  rm(full)
  
  #BINWMEANCP.DF[[chr]] <- cbind.data.frame(chr=chr, bin=1:totbin, wmeanCp0, wmeanCp)
  BINWMEANCP.DF <- cbind.data.frame(bin=1:totbin, wmeanCp0, wmeanCp, stringsAsFactors=FALSE)
  rownames(BINWMEANCP.DF) <- NULL
  save(BINWMEANCP.DF, file=paste0(out.dir, "/", gcb, "_", chr, "_weightedMeanCp.RData"))
  
  #rm(totbin)
  rm(totbin, BINWMEANCP.DF); gc()
  
}

#BINWMEANCP.DF <- do.call("rbind.data.frame", c(BINWMEANCP.DF, stringsAsFactors=FALSE))
#rownames(BINWMEANCP.DF) <- NULL
#save(BINWMEANCP.DF, file=paste0(out.dir, "/", gcb, "_weightedMeanCp.RData"))

end.time <- Sys.time()
end.time-start.time 

# rm(list=ls()); gc()




