################################################################################
# Generate per chr the IJ.RT.TYPES object for contact-wise plot of RT vs. Cp.
# IJ.RT.TYPES contains consensus value per contact depending on ij.fnx(na.rm=F)
# for all rt types i.e. all, nontumor, tumor (so 3 columns). 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/20_ReplicationTiming")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
src.dir = paste0(wk.dir, "/out_plotdata")
out.dir = paste0(wk.dir, "/out_contact_Cp_plotdata")
### OTHER SETTINGS #############################################################
rt.calc = "mean" # CALCREPLACE
rt.types = c("all", "nontumor", "tumor")
chrs = paste0("chr", c(21:22))
persist.id = "Persist_min2Mb"
src.id = paste0(rt.calc, "_normNot1_setpointcountGrEq3_40000bpHiCres")
ij.fnx = "mean" # "median" # IJFNXREPLACE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(chr in chrs){
  
  # Load contact data
  
  load(paste0(persist.dir, "/", chr, "_", persist.id, ".RData"))
  ij.mx <- cbind(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, Cp=PERSIST.MX$ntis)
  rm(PERSIST.MX)
  rownames(ij.mx) <- NULL
  
  vals.mx <- ij.mx[,1:2]
  
  IJ.RT.TYPES <- foreach(rt.type=rt.types, .inorder=T, .combine="cbind") %do% {
      
    # Load source data (rt plot data)
    load(paste0(src.dir, "/", chr, "_", rt.type, "_", src.id, ".RData"))
    
    # Calculate contact value by taking ij.fnx(region1,region2)
    
    vals.mx[,1] <- vals[as.integer(ij.mx[,1])]
    vals.mx[,2] <- vals[as.integer(ij.mx[,2])]
    
    if(ij.fnx == "mean"){
      ij.vals <- rowMeans(vals.mx, na.rm=F)
    } else {
      
      eval(parse(text=paste0(
        'ij.vals <- apply(vals.mx, MARGIN=1, FUN=', ij.fnx, ', na.rm=F)'
      )))
      
    }
    
    return(ij.vals)
    
  } # rt.types foreach loop end
  
  dimnames(IJ.RT.TYPES)[[2]] <- rt.types
  IJ.RT.TYPES <- cbind(IJ.RT.TYPES, Cp=ij.mx[,"Cp"])
  
  save(IJ.RT.TYPES, file=paste0(out.dir, "/", chr, "_ij.fnx", ij.fnx, "_allrttypes_", src.id, ".RData"))
  
} # chrs for loop end

# rm(list=ls()); gc()