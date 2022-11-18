################################################################################
# Combine CII.MX of all chromosomes
# deva, R/3.6.0-newgcc, gcc/4.9.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity")
lib = paste0(home.dir, "/DPhil/lib")
compl.dir = paste0(wk.dir, "/out_constraints_GfreeSingleNorm/merged_final")
out.dir = compl.dir
### OTHER SETTINGS #############################################################
chr.v = "chr18" #paste("chr", c("ALL", 1:22, "X"), sep="")
combineChr = FALSE
gcb = "min2Mb"
kmer.len = 7L
bin.len = 4e4L
type = "kmer"
affix = ""
gap.rng = c(50,100) # j - i - 1, closed range, Set to NULL if not filtering
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(RColorBrewer)
source(paste0(lib, "/HiCHybridPlot.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( !is.null(gap.rng) ){
  out.dir <- paste0(out.dir, "/filter_gap_jMINUSiMINUS1equals", paste(gap.rng, collapse="To"), "bins")
  if( !dir.exists(out.dir) ){ dir.create(out.dir) }
}

if(type=="align"){
  ylab <- list(`C||`=bquote(bold( "c"["||"]~"align" )) )
} else if(type=="kmer"){
  ylab <- list(`C||`=bquote(bold( "c"["||"]~"kmer" )),
               Gfree=bquote(bold( "Gfree, kcal/mol" )),
               sdDifference=bquote(bold( "s ("~"c"["||"]~")" ))
  )
} else {
  stop("Invalid type.")
}
#---------------------------------------
if(combineChr){
  CII.MX <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
    # Load CII.MX
    load(file=paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"))
    # Filter for given gcb
    CII.MX <- CII.MX[!is.na(CII.MX[,"Cp"]), ]
    print(paste0(chr, " data obtained!"), quote=FALSE)
    return(CII.MX)
  })
  CII.MX <- do.call("rbind", CII.MX)
  save(CII.MX, file=paste0(out.dir, "/chrALL_", type, "_", gcb, affix, ".RData"))
  chr.v <- "chrALL"
  print("Combined chr data.", quote=FALSE)
  
} 
#---------------------------------------
for(chr in chr.v){
  if(!combineChr){
    load(file=paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"))
  }
  
  if( !is.null(gap.rng) ){
    print(paste0(out.dir, ": Filtering contacts based on given closed gap range."))
    gaps <- CII.MX[,"j"] - CII.MX[, "i"] - 1
    CII.MX <- CII.MX[ gaps >= gap.rng[[1]] & gaps <= gap.rng[[2]] & 
                      !is.na(CII.MX[,"Cp"]), ]
  }
  
  cp.v <- unname(CII.MX[,"Cp"])
  colnme <- colnames(CII.MX)[!colnames(CII.MX)%in%c("i", "j", "Cp", "group")]
  CII.MX <- CII.MX[,colnme, drop=F]
  dimnames(CII.MX) <- list(NULL, colnme)
                   
  ylab <- ylab[colnames(CII.MX)]
  HiCHybridPlot(
    out.dir=out.dir,
    out.name=paste0(chr, "_", type, "_", gcb, affix),
    HYB.MX=CII.MX,
    # Corresponding order to HYB.MX (column)
    cp.v=cp.v,
    # Corresponding order to HYB.MX (row)
    label.list=ylab
  )
  rm(CII.MX, cp.v); gc()
  print(paste0(chr, " done!"), quote=FALSE)
}

# rm(list=ls())






