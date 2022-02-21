################################################################################
# Calculate p-value for c|| vs Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/11_Complementarity")
    CII.dir  = paste0(wk.dir, "/z_ignore_git/out_constraints/merged_final")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" 
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/11_Constraints")
    CII.dir  = paste0(wk.dir, "/out_constraints/merged_final")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")

out.dir  = paste0(wk.dir, "/out_calc_significance")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chrALL"
type.v = c("Gfree", "sdDifference") #c("kmer", "align") #, "Gfree", "sdDifference")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
source(paste0(lib, "/doCorTest.R"))
source(paste0(lib, "/doVarTest.R"))
source(paste0(lib, "/compareManyDist.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(type in type.v){
  
  if( type %in% c("Gfree", "sdDifference") ){
    
    dtype.id <- "kmer"
    val.id <- type
    if(type=="Gfree"){ mult <- 10^5 } else { mult <- 10^4 }
    
  } else if( type %in% c("kmer", "align") ){
    
    dtype.id <- type
    val.id <- "C||"
    mult <- 1
    
  } else {
    stop("Invalid argument.")
  } 
  
  fle.id <- paste0(chr, "_", dtype.id, "_", gcb)
  load(paste0(CII.dir, "/", fle.id, ".RData"))
  
  CII.MX <- CII.MX[! (is.na(CII.MX[,val.id]) | is.na(CII.MX[,"Cp"])), ]
  
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
  #doCorTest(xval=Cps, yval=vals, alt="two.sided", exactpval=F, out.dir, out.name)
  
  Cps <- factor(x=as.character(Cps), 
                levels=as.character(sort(unique(Cps)))
                )
  
  # ANOVA, Kruskal-Wallis H-test
  doVarTest(xval=vals, grp=Cps, out.dir, out.name)
    
  ## Post-hoc: Pairwise t-test and MWw test; takes longer
  #compareManyDist(xval=vals, grp=Cps, alt="two.sided", out.dir, out.name)  
    
  print(paste0(out.name, " done!"), quote=F)
  
  rm(vals, Cps, type, dtype.id, val.id, fle.id, out.name)
  gc()
  
}

# rm(list=ls()); gc()