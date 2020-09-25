################################################################################
# Bin data
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/persist_HiCNorm")
out.dir = paste0(wk.dir, "/binmx/out_bindata_1perc_HiCNorm")
### OTHER SETTINGS #############################################################
gcb = "min2Mb" # "min2Mb" | "min05Mb"
chr.v = paste0("chr", c(1:22, "X"), sep="")
nCPU = 4L
CT.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
         "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
Cp.v = 1:21
Cs.cutoff = 0.01 # Target fraction of contacts with the highest Cs
bin.len = 40000
# If matchCp21N=TRUE, override Cs.cutoff and the code will match number of 
# contacts with highest Cs to number of Cp=max(Cp.v) contacts per chr.
# Cs.cutoff marker will be 2.
matchMaxCpN = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, "..."), quote=FALSE)
# Chromosome length file
chrLen.df <- read.table(file=chrLenfile, stringsAsFactors=FALSE, header=TRUE,
                        colClasses=c("character", "integer", "integer"))
Cp.v <- sort(unique(Cp.v), decreasing=FALSE)
CT.v <- sort(unique(CT.v))
col.nme <- paste("s_Cp_", rep(Cp.v, each=length(CT.v)), 
                 "_ct_", rep(CT.v, times=length(Cp.v)), 
                 "_e", sep="")
if( any(duplicated(col.nme)) ){ stop("Checkpoint 1.") }
col.nme <- c("start", "end", col.nme)

toExport <- c("persist.dir", "out.dir", "gcb", "chr.v", "CT.v", "Cp.v", 
              "Cs.cutoff", "bin.len", "chrLen.df", "col.nme", "matchMaxCpN")
chr.v.len <- length(chr.v)
#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
        .inorder=FALSE, .export=toExport, 
        .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  for(i in itr){
    
    chr <- chr.v[i]
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    PERSIST.MX$valsum <- PERSIST.MX$control <- NULL
    chr.len <- chrLen.df[chrLen.df$chromosome==chr,"length.bp"]
    bin.last <- ceiling(chr.len/bin.len)
    
    # Initialise matrix
    BIN.MX <- matrix(data=0, nrow=bin.last, ncol=length(col.nme), 
                     dimnames=list(1:bin.last, col.nme))
    BIN.MX[,"end"] <- (1:bin.last)*bin.len
    BIN.MX[,"start"] <- (BIN.MX[,"end"])-bin.len+1
    BIN.MX[bin.last,"end"] <- chr.len; rm(chr.len, bin.last)
    
    for(ct in CT.v){
      ct.TF <- PERSIST.MX$hits[[ct]]>0
      ij.df <- cbind(PERSIST.MX$hits[ct.TF,c("i", "j", ct)],
                     Cp=PERSIST.MX$ntis[ct.TF])
      rm(ct.TF)
      ij.df <- ij.df[order(ij.df[[ct]], decreasing=TRUE),]
      if(matchMaxCpN){
        ij21.len <- sum(ij.df$Cp==max(1:21)) 
        ij.df[1:ij21.len,ct] <- 2; rm(ij21.len)
        Cs.cutoff <- 2
      } else {
        # Mark top contacts based on Cs
        toplen <- ceiling(nrow(ij.df)*Cs.cutoff)
        ij.df[ ij.df[[ct]] > ij.df[toplen,ct], ct] <- Cs.cutoff
        count <- sum(ij.df[[ct]]==Cs.cutoff)
        print(paste0(chr, "_", ct, " Target: ", toplen, " Actual: ", count, 
                     " Difference: ", toplen-count), quote=FALSE)
        rm(toplen)
      }
      
      for(Cp in Cp.v){
        Cp.TF <- ij.df$Cp==Cp
        binsCTCpCs <- unique(unlist(
          ij.df[ ij.df[[ct]]==Cs.cutoff & Cp.TF,c("i","j") ]
        ))
        binsCTCp <- unique(unlist( ij.df[Cp.TF,c("i", "j")] ))
        rm(Cp.TF)
        nme <- paste0("s_Cp_", Cp, "_ct_", ct, "_e")
        BIN.MX[binsCTCpCs,nme] <- Cs.cutoff
        BIN.MX[setdiff(binsCTCp, binsCTCpCs),nme] <- 1
        rm(nme); gc()
      }
      rm(ij.df); gc()
      print(paste0(ct, " done!"), quote=FALSE)
    }  # CT.v for loop end
    
    rm(PERSIST.MX); gc()
    save(BIN.MX, file=paste0(out.dir, "/", chr, "_", gcb, "_bindata.RData"))
    rm(BIN.MX); gc()
    print(paste0(chr, " done!"), quote=FALSE)
    
  } # itr for loop end
}

# rm(list=ls()); gc()



