################################################################################
# Compare set of contacts present in a tissue with set without tissue filtering.
# Define contacts to be compared based on contact gap. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
#persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
persist.dir = paste0(wk.dir, "/out_basePersist")
### OTHER SETTINGS #############################################################
gcb = "min2Mb" # "min2Mb" | "min05Mb"
chr.v = "chr1" #paste0("chr", c(1:22, "X"), sep="")

gap.bin = 650
cp.v = 19:21
ct.v = "FC" #c("FC", "ESC", "LC")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
cp.v <- unique(cp.v)
print(paste0("Cp: ", paste(cp.v, collapse=",")), quote=F)

for(chr in chr.v){
  
  for(ct in ct.v){
    
    #load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, "_topCP4.RData"))
    
    out.id <- paste0(gcb, "_", chr, "_gapBin", gap.bin, "_", ct)
    if( !ct%in%colnames(PERSIST.MX$hits) ){
      stop( paste0(out.id, ": ", ct, " not in PERSIST.MX") )
    }
    
    if( any(!cp.v%in%PERSIST.MX$ntis) ){
      stop( paste0(out.id, ": Not all elements of cp.v in PERSIST.MX") )
    }
    
    gap.TF <- (PERSIST.MX$hits$j-PERSIST.MX$hits$i) > gap.bin
    cp.TF <- PERSIST.MX$ntis%in%cp.v
  
    all.TF <- cp.TF & gap.TF
    ct.TF <- (PERSIST.MX$hits[[ct]] > 0) & all.TF
    
    if( identical(all.TF, ct.TF) ){
      
      print(paste0(out.id, " contact set identical to set not cell/tissue-filtered."),
            quote=F)
    
    } else {
      
      print(paste0(out.id, " contact set different from set not cell/tissue-filtered."),
            quote=F)
      
    }
   
  } # ct.v for loop end
  
} # chr.v for loop end
         
# rm(list=ls()); gc()