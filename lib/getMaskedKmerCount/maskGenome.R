################################################################################
# Function for masking one or any combination of genomic features; IRanges based
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(IRanges)
################################################################################
################################################################################
maskGenome <- function(chr = "chr21",
                       seq="character vector or string", 
                       maskbed = "bed-format table", 
                       maskingChar = "m",
                       reduce = FALSE,
                       split = FALSE
                       ){

  # Subset data for the right chromosome
  maskbed <- maskbed[maskbed[,1]==chr,]

  # Simplify ranges
  if(reduce==TRUE){
    
    maskbed <- as.data.frame(reduce(
      IRanges( start=as.numeric(maskbed[,2]),
               end=as.numeric(maskbed[,3]) )
    )) 
    # Remove width column
    maskbed <- cbind(chr, maskbed[,-3])
    
  }
  
  # Mask with maskingChar
  len <- length(maskbed[,1])
  
  for(i in 1:len){
    
    if(split){
      seq[maskbed[i,2]:maskbed[i,3]] <- maskingChar
    } else{
      m.len <- maskbed[i,3]-maskbed[i,2] + 1L
      substr(x=seq, 
             start=maskbed[i,2], 
             stop=maskbed[i,3]) <- paste(rep(x=maskingChar, times=m.len), collapse="")
    }
    
  } # len for loop end
  
  return(seq)
  
}

# rm(list=ls())
