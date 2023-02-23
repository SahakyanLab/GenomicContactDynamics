################################################################################
# Title
# deva, R/3.5.0-newgcc, gcc/4.9.2
# deva, R/3.6.0-newgcc, gcc/4.9.2
# Mac, R/3.5.2, R/3.6.1
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
checkANDfilterRTdata <- function(rt, HiC.res, chrs, chrLen.file){
  
  chrLen.df <- read.delim(file=chrLen.file, header=T)
  
  for( chr in unique(rt$chroms) ){
    
    # Check that rt data resolution equal to Hi-C resolution
    
    rt.res <- unique(diff( rt$starts[rt$chroms==chr] ))
    if( !identical(rt.res, HiC.res) ){
      warning( paste0("checkANDfilterRTdata(): ", chr, ": rt.res and HiC.res not identical.") )
    }
    
    # Check if rt data has all hg19 bins and in order 
   
    rt$bin <- ceiling(rt$starts / rt.res) 
    chr.len <- chrLen.df$length.bp[chrLen.df$chromosome==chr]
    tot.bin <- ceiling(chr.len / rt.res)
    
    chr.TF <- rt$chroms==chr
    
    if(
      !identical(as.numeric(rt$bin[chr.TF]), as.numeric(1:tot.bin))
    ){
      warning( paste0("checkANDfilterRTdata(): ", chr, ": Missing chr bins in rt.") )
    }
    
    if( max(rt$starts[chr.TF]) > chrLen.df$length.bp[chrLen.df$chromosome==chr] ){
      warning( paste0(chr, ": max rt start coordinate > chr length.") )
    }
    
  }
  
  # # Filter, set to NA those bins not satisfying criteria
  # 
  # col.TF <- !colnames(rt)%in%c("chroms", "starts")
  # dropBin.TF <- rt$set.count < 3 | rt$point.count < 3 | rt$norm == 1 
  # rt[dropBin.TF,col.TF] <- NA
  # 
  # message("checkANDfilterRTdata(): rt data checked and filtered.")
  
  message("checkANDfilterRTdata(): rt data checked.")
  
  return(rt)

}

# rm(list=ls()); gc()