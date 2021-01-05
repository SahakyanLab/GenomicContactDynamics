################################################################################
# Make dataframe of Cs, Cp, CII (kmer and align) of all contacts per cell/tissue
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
melt.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/persist_HiCNorm")
compl.dir = paste0(wk.dir, "/out_constraints/merged_final")
#compl.dir = paste0(wk.dir, "/out_constraints")
out.dir = paste0(wk.dir, "/out_compare_HiCNormCs")
### OTHER SETTINGS #############################################################
ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
         "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "LC", "FC")
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
nCPU = 2L #~15G
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
for(ct in ct.v){
  id <- paste0("chrALL_", gcb, "_", ct)
  chr.v.len <- length(chr.v)
  toExport <- c("persist.dir", "compl.dir", "ct", "chr.v", "gcb")
  
  #### PARALLEL EXECUTION #########
  CSCPCII.MX <- foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
                        .inorder=TRUE, .combine="rbind",
                        .export=toExport, .noexport=ls()[!ls()%in%toExport]
  ) %op% {
    chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
      chr <- chr.v[i]
      #-------------------
      # Load MELT.MX
      load(file=paste0(melt.dir, "/human_", chr, "_allcontacts.RData"))
      MELT.MX$upper.tri.nocontact <- MELT.MX$rest <- NULL
      ij.ct.TF <- MELT.MX$upper.tri[,ct]!=0
      # Initialize output
      df <- data.frame(chr=chr,
                       Cs=MELT.MX$upper.tri[ij.ct.TF,ct],
                       Cp=NA_integer_,
                       CIIkmer=NA,
                       CIIalign=NA)
      rownames(df) <- rownames(MELT.MX$upper.tri)[ij.ct.TF]
      rm(MELT.MX, ij.ct.TF); gc()
      #-------------------
      # Load PERSIST.MX 
      load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
      ij.cpct.TF <- PERSIST.MX$hits[,ct]!=0
      ij.cp.rwnme <- rownames(PERSIST.MX$hits[ij.cpct.TF,])
      df[ij.cp.rwnme, "Cp"] <- PERSIST.MX$ntis[ij.cpct.TF]
      rm(PERSIST.MX, ij.cp.rwnme); gc()
      #-------------------
      # Load CII.MX kmer
      load(file=paste0(compl.dir, "/", chr, "_kmer_", gcb, ".RData"))
      df$CIIkmer <- CII.MX[rownames(df),"C||"]
      rm(CII.MX); gc()
      #-------------------
      # Load CII.MX align
      load(file=paste0(compl.dir, "/", chr, "_align_", gcb, ".RData"))
      df$CIIalign <- CII.MX[rownames(df),"C||"]
      rm(CII.MX); gc()
      #-------------------
      print(chr, quote=FALSE)
      return(df)
    })
    return( do.call("rbind", chunk) ) 
  }
  ### END OF PARALLEL EXECUTION ###
  save(CSCPCII.MX, file=paste0(out.dir, "/", id, "_CsCpCII.RData"))
  rm(CSCPCII.MX); gc()
  print(paste0(ct, " done!"), quote=FALSE)
}

# rm(list=ls())
