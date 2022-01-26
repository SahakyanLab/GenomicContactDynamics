################################################################################
# Extract each matrix and save as document type with .mat extension
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(data.table)
### FUNCTION ###################################################################
extractEachMatrix <- function(
  
  matrix.filepath = 'source matrix',
  chrlen.filepath = 'output file path',
  out.dir = 'output dir', 
  out.prefix = '',
  out.suffix = '',
  bin.len = 'bin resolution',
  chr.ids = '1:12'
  
){
  
  source.mat <- data.table::fread(file=matrix.filepath, header=F, stringsAsFactors=F, 
                                  data.table=F, sep="\t")
  mat <- data.matrix(source.mat)
  rm(source.mat)
  gc()
  dimnames(mat) <- NULL
  
  tot.bin <- unique(dim(mat))
  
  if( !isSymmetric(mat) ){
    stop("mat not symmetric.")
  }
  
  chrlen.df <- read.delim(file=chrlen.filepath, header=T, stringsAsFactors=F)
  
  chr.ids <- as.character(chr.ids) 
  chr.bins <- setNames(obj=rep(NA_integer_, times=length(chr.ids)),
                       nm=chr.ids)
  
  start.bin <- 1
  for(chr.id in chr.ids){
    
    chr.bins[chr.id] <- ceiling( chrlen.df$length.bp[chrlen.df$chromosome==paste0("chr", chr.id)] / bin.len )
    print(paste0("chr ", chr.id, ": ", chr.bins[[chr.id]], " expected!"), quote=F)
    
    end.bin <- start.bin + chr.bins[[chr.id]] - 1
    print(paste0("chr ", chr.id, " coordinates: ", start.bin, ":", end.bin), quote=F)
    mat.chr <- mat[ start.bin:end.bin, start.bin:end.bin ]
    
    # start.bin for next chr.id
    start.bin <- end.bin + 1
    
    mat.chr.bin <- unique(dim(mat.chr))
    if( mat.chr.bin != chr.bins[[chr.id]] ){
      print(paste0("chr ", chr.id, ": Wrong matrix dimension!"), quote=F)
    } else {
      
      chr.bins[chr.id] <- mat.chr.bin
      
      write.table(x=mat.chr, 
                  file=paste0(out.dir, "/", out.prefix, chr.id, out.suffix, ".mat"), 
                  sep="\t", quote=F, row.names=F, col.names=F)
      
    }
    
    if( !isSymmetric(mat.chr) ){
      stop(paste0("chr ", chr.id, ": mat.chr not symmetric."))
    } 
    
    rm(mat.chr)
    gc()
    
    print(paste0("chr ", chr.id, " done!"), quote=F)
    
  }
  
  if(tot.bin != sum(chr.bins) ){
    stop("Number of bins does not sum up.")
  }
  
}

# rm(list=ls()); gc()