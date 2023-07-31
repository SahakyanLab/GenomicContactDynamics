################################################################################
# Extract dense matrix from .hic (currently) containing sparse 3-column contact
# data. 
# A normalised data from straw has low read contacts removed by setting to NaN.
# When there is some form of normalisation method specified (i.e. when norm.method
# is not "NONE"), the code expects this and sets NaN to 0, which can be set to NA 
# when plotting. Then it checks again for non-finite values and then displays 
# final count range. Final output should have no non-finite values.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(strawr)
### FUNCTION ###################################################################
extractDenseMxMap <- function(chr,
                              bin.len,
                              norm.method = "KR",
                              mx.type = "oe",
                              # "symmetric", "upper", "lower" - all with diagonal values
                              mx.structure = "symmetric", 
                              straw.unit = "BP",
                              hic.file,
                              chrlen.file,
                              out.dir,
                              out.name
                              ){
  
  message(paste0("extractDenseMxMap(): Processing chr", chr, "..."))
  
  chrlen.df <- read.delim(file=chrlen.file, sep="\t", stringsAsFactors=F, header=T)
  
  # 
  hic.bpres.avail <- readHicBpResolutions(fname=hic.file)
  
  if(! bin.len %in% hic.bpres.avail){
    stop("extractDenseMxMap(): Required bin.len not available in .hic")
  }
  
  # 
  hic.chrs.avail <- readHicChroms(fname=hic.file)$name
  if(!chr %in% hic.chrs.avail){
    stop("extractDenseMxMap(): Not all chrs in hic available chrs.")
  }
  
  # bin 1 start - bin 2 start - cf
  melted.df <- straw(norm=norm.method, fname=hic.file, chr1loc=chr, chr2loc=chr, 
                     unit=straw.unit, binsize=bin.len, matrix=mx.type)
  melted.mx <- data.matrix(melted.df)
  
  # Check values in melted.mx
  
  chr.len <- chrlen.df$length.bp[chrlen.df$chromosome == paste0("chr", chr)]
  if( max(melted.mx[,1:2]) >= chr.len ){
    rm(melted.mx)
    stop("extractDenseMxMap(): Max start coordinate in melted.mx exceed or equal chr length.")
  }
  
  if( min(melted.mx[,3], na.rm=T) < 0 ){
    stop(paste0("extractDenseMxMap(): Negative counts in .hic data."))
  }
  
  # Downstream code cannot deal with non-finite count. Normalisation (KR) removes low read counts by 
  # setting to NaN. Assign them as 0, which will be set to NA when plotting. Other mx.type should 
  # have all finite counts.
  if(norm.method != "NONE"){
    
    nan.TF <- is.nan(melted.mx[,3])
    melted.mx[,3][nan.TF] <- 0
    
    message(paste0("extractDenseMxMap(): ", sum(nan.TF), " NaNs set to 0."))
    
  }
  
  if( any(!is.finite(melted.mx[,3])) ){
    stop(paste0("extractDenseMxMap(): Non-finite counts in melted.mx even after removing NaNs from normalisation."))
  } 
  
  message(paste0("extractDenseMxMap(): Final count range is ", paste(range(as.numeric(melted.mx[,3]), na.rm=F), collapse="-")) )
  
  # straw mx 0-based so convert to 1-based as expected by downstream code
  
  gap.bp <- melted.mx[,2] - melted.mx[,1]
  if( all(gap.bp >= 0) ){
    melted.mx[,1] <- melted.mx[,1] + 1
    melted.mx[,2] <- melted.mx[,2] + 1
  } else{
    stop("extractDenseMxMap(): Negative contact gaps in melted.mx.")
  }
  
  # straw returns start coordinates in bp so convert to bin number
  
  melted.mx[,1] = ceiling(melted.mx[,1] / bin.len)
  melted.mx[,2] = ceiling(melted.mx[,2] / bin.len)
  
  # straw() returns sparse 3-column contact data frame so convert to contact dense matrix
  
  chr.bin.num <- ceiling(chr.len / bin.len)
  dense.mx <- matrix(0, nrow=chr.bin.num, ncol=chr.bin.num)
  
  # https://stackoverflow.com/questions/59536834/is-it-faster-to-fill-a-matrix-by-row-or-to-transpose-a-matrix-filled-by-columns
  if(mx.structure == "lower"){
    dense.mx[melted.mx[,c(2,1)]] = melted.mx[,3]
  } else if(mx.structure == "upper"){
    dense.mx[melted.mx[,c(1,2)]] = melted.mx[,3]
  } else if(mx.structure == "symmetric"){
    
    dense.mx[melted.mx[,c(2,1)]] = melted.mx[,3] # lower triangle
    dense.mx[melted.mx[,c(1,2)]] = melted.mx[,3] # upper triangle
    
    if( !isSymmetric(dense.mx) ){
      stop("extractDenseMxMap(): dense.mx not symmetric.")
    }
 
  } else {
    stop("extractDenseMxMap(): Invalid mx.structure argument.")
  }
   
  ## Make symmetric 
  #dense.mx[lower.tri(dense.mx)] <- t(dense.mx)[lower.tri(dense.mx)]
  
  out.file <- paste0(out.dir, "/", out.name)
  write.table(x=dense.mx, file=out.file, sep="\t", quote=F, row.names=F, col.names=F, na="NA")
  
}

# rm(list=ls()); gc()