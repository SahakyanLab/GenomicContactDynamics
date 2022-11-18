################################################################################
# Calculate p-value for c|| vs Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
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
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity")
constraints.id = c("GfreeSingleNorm", "hg19_rm_GfreeSingleNorm")[[arr1.repl]]
CII.dir = paste0(wk.dir, "/out_constraints_", constraints.id, "/merged_final")
out.dir  = paste0(wk.dir, "/out_calc_significance_", constraints.id)
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr22"
type.v = c("kmer", "Gfree", "sdDifference", "align")
mult.v = setNames(c(1, 1, 10^4, 1), nm=type.v) # Multiplier before doing test, for less memory usage
affix = ""
gap.rng = list(c(50,50), c(50,100))[[arr2.repl]] # j - i - 1, closed range, Set to NULL if not filtering
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
source(paste0(lib, "/doCorTest.R"))
source(paste0(lib, "/doVarTest.R"))
source(paste0(lib, "/compareManyDist.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( !is.null(gap.rng) ){
  out.dir <- paste0(out.dir, "/filter_gap_jMINUSiMINUS1equals", 
                    paste(gap.rng, collapse="To"), "bins")
  if( !dir.exists(out.dir) ){ dir.create(out.dir) }
}

for(type in type.v){
  
  mult <- mult.v[[type]]
  if( type %in% c("Gfree", "sdDifference") ){
    
    dtype.id <- "kmer"
    val.id <- type
    #if(type=="Gfree"){ mult <- 10^5 } else { mult <- 10^4 }
    
  } else if( type %in% c("kmer", "align") ){
    
    dtype.id <- type
    val.id <- "C||"
    #mult <- 1
    
  } else {
    stop("Invalid argument.")
  } 
  
  fle.id <- paste0(chr, "_", dtype.id, "_", gcb, affix)
  try(load(paste0(CII.dir, "/", fle.id, ".RData")))
  
  is_valid_ij <- ! (is.na(CII.MX[,val.id]) | is.na(CII.MX[,"Cp"]))
  if( !is.null(gap.rng) ){
    gaps <- CII.MX[,"j"] - CII.MX[, "i"] - 1
    is_valid_ij <- is_valid_ij & (gaps >= gap.rng[[1]] & gaps <= gap.rng[[2]])
    rm(gaps)
  }
  CII.MX <- CII.MX[is_valid_ij,]
  
  vals <- CII.MX[,val.id] * mult
  if( type %in% c("Gfree", "sdDifference") ){
    
    print(range(vals, na.rm=T), quote=F)
    vals <- round(vals, digits=4)
    print(range(vals, na.rm=T), quote=F)
    
  }
    
  Cps <- CII.MX[,"Cp"]
  
  rm(CII.MX)
  
  out.name <- paste0(fle.id, "_", val.id) 
  
  # Correlation 
  doCorTest(xval=Cps, yval=vals, alt="two.sided", exactpval=F, out.dir, out.name)
  
  Cps <- factor(x=as.character(Cps), 
                levels=as.character(sort(unique(Cps)))
                )
  
  # ANOVA, Kruskal-Wallis H-test
  doVarTest(xval=vals, grp=Cps, out.dir, out.name)
    
  # Post-hoc: Pairwise t-test and MWw test; takes longer
  compareManyDist(xval=vals, grp=Cps, alt="two.sided", out.dir, out.name)  
    
  print(paste0(out.name, " done!"), quote=F)
  
  rm(vals, Cps, type, dtype.id, val.id, fle.id, out.name)
  gc()
  
}

# rm(list=ls()); gc()