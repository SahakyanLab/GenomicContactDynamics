################################################################################
# Function for masking a few types of genomic features and any of their combination
# bedtools based
# REQUIRES:
# library(bedr) # to access bedtools via R
# library(data.table)
# library(foreach)
# masks: repeat | CTCFBS | CpG | transcript
# mask file should be in bed file format
# mask file format: genomeVer_mask e.g. hg19_repeat
################################################################################
################################################################################
maskGenome <- function(genomeObj=obj, genomeVer="hg19", chr="chr21", 
                       maskType=c("repeat", "CpG"), 
                       maskingLib="/Users/ltamon/Database/maskingLib",
                       maskingChar="n",
                       bedtoolsPATH="/Users/ltamon/prog/bedtools2/bin"
                       ){
  
  library(bedr) # to access bedtools via R
  library(data.table)
  library(foreach)
  
  # add path containing bedtools to PATH (temporary)
  PATHorig <- Sys.getenv("PATH")
  if( !grepl(pattern=bedtoolsPATH, x=Sys.getenv("PATH")) ){
    Sys.setenv(PATH=paste(PATHorig, bedtoolsPATH, sep=":"))
  }
  # confirm that bedtools can be accessed
  if(!check.binary("bedtools", verbose=TRUE)){
    stop("bedtools not found.")
  }
  
  # merge masks
  mask.all <- foreach(mask=maskType, .inorder=TRUE, .combine="c"
  ) %do% {
    maskfile=paste0(maskingLib, "/", dir(maskingLib, pattern=mask))
    maskbed <- fread(file=maskfile, data.table=FALSE, stringsAsFactors=FALSE, 
                    select=c("chr", "start", "end"))
    # subset by chr
    maskbed <- maskbed[maskbed$chr==chr,]
    bed2index(maskbed)
  }
  mask.all <- bedr.merge.region( bedr.sort.region(mask.all) )
  rm("maskfile", "maskbed"); gc()
  
  # convert mask to data frame
  mask.all <- index2bed(mask.all)
  
  # mask with maskingChar
  len <- length(mask.all[,1])
  for(i in len){
    genomeObj$seq[mask.all[i,1]:mask.all[i,1]] <- maskingChar
  }
  rm("mask.all", "len"); gc()

  # revert PATH variable to original state
  Sys.setenv(PATH=PATHorig)
  
  return(genomeObj)
  
}

# rm(list=ls())
