################################################################################
# Extract each fasta from a multiple-fasta file. Use make.uppercase argument to 
# change all bases to upper case in final fasta. After changing (or leaving case
# as is), change invalid bases not in acceptable.bases to N (this is a 
# case-sensitive action done AFTER applying make.uppercase argument). Regardless
# of case in the final fasta loadGenome() was made to return upper cases by default
# but can be modified through the case argument of readfasta() used by loadGenome(). 
# Create file with chromosome lengths.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# lib = paste0("/Users/ltamon/DPhil/lib")
# library(data.table)
# source(paste0(lib, "/TrantoRextr/GEN_readfasta.R"))                              
# source(paste0(lib, "/TrantoRextr/UTIL_readLinesFast.R"))  
# source(paste0(lib, "/TrantoRextr/GEN_loadGenome.R"))
### FUNCTION ###################################################################
extractEachFasta <- function(
  
  fasta.filepath = 'source file of fastas',
  out.dir = 'output dir',
  chr.ids = 'chromosomes in order based on source fasta, e.g. 1:2',
  make.uppercase = 'convert bases to upper case in final fasta, if FALSE, 
  # keep case as is',
  acceptable.bases = "case-aware; change to N the bases not acceptable; filtering
  # done after modifying case (if applicable)",
  
  fastafile.ending = 'fasta file extension',
  genome.prefix = 'genome.prefix of individual fasta',
  
  # Create file with chromosome lengths and number of bins given bin.length
  chrlen.filepath = 'output file path',
  bin.len = 50000
  
){

  source.fas <- readLines(con=fasta.filepath)
  
  # Check for and remove blank lines if present
  is.blank <- !grepl(x=source.fas, pattern="[[:alnum:]]", fixed=F)
  if( any(is.blank) ){
    
    source.fas <- source.fas[!is.blank] 
    warning("Blank lines in fasta file removed!")
    
  }
  rm(is.blank)
  
  header.ind <- grep(x=source.fas, pattern=">", fixed=T) 
  seq.start.ind <- header.ind + 1
  seq.end.ind <- c( header.ind[-1] - 1, length(source.fas) )
  seq.len <- length(header.ind)
  
  if( length(chr.ids) != length(header.ind) ){
    stop("Length of chr.ids argument not equal to number of fasta.")
  }
  
  acceptable.bases.pattern <- paste(acceptable.bases, collapse="|")
  acceptable.bases.pattern <- paste0("[^", acceptable.bases.pattern, "]")
    
  # Save the individual fasta
  for(i in 1:seq.len){
    
    fas.seq <- source.fas[ seq.start.ind[i]:seq.end.ind[i] ]
    
    if(make.uppercase){
      
      print("Converting to upper case.")
      fas.seq <- toupper(fas.seq)
      
    }
      
    fas.seq <- gsub(x=fas.seq, pattern=acceptable.bases.pattern, 
                    replacement="N", fixed=F, ignore.case=F)
    
    fas.seq <- c(source.fas[header.ind[i]], fas.seq)
    
    out.name <- paste0(genome.prefix, chr.ids[i], fastafile.ending)
    writeLines(fas.seq, con=paste0(out.dir, "/", out.name), sep="\n")
    
    print(paste0(fas.seq[1], " -> ", out.name, " extracted!"), quote=F)
    
    rm(fas.seq, out.name)
    gc()
    
  }
  
  # 
  chr.lens <- setNames(object=rep(x=NA, times=length(chr.ids)), nm=chr.ids)
  
  for(chr.id in chr.ids){
    
    # Set to return uppercase letters
    loadGenome(PATH.genome=out.dir, genome.prefix=genome.prefix,
               fastafile.ending=fastafile.ending, chr.id=chr.id, # any of c(1:5,"MT", "chloroplast")
               silent=F, remove.other.loads=F, split=T)
    
    eval(parse(text=paste0(
      'genome.obj <- ', genome.prefix, chr.id, fastafile.ending
      #, "; rm(",genome.prefix, chr.id, fastafile.ending, ")"
    )))
    
    chr.lens[chr.id] <- genome.obj$length
    
    print(paste0(chr.id, " bases tally: "), quote=F)
    print(table(genome.obj$seq, useNA="ifany"), quote=F)
    
    if( any( ! genome.obj$seq %in% unique(toupper(acceptable.bases)) ) ){
      stop(paste0(chr.id, ": Invalid base."))
    } 
    
    rm(genome.obj)
    gc()
    
    print(paste0(chr.id, " checked!"), quote=F)
    
  }
  
  #
  chr.lens <- stack(chr.lens)[,2:1]
  colnames(chr.lens) <- c("chromosome", "length.bp")
  chr.lens$chromosome <- paste0("chr", chr.lens$chromosome)
  chr.lens[[paste0("bins.", bin.len, "bp")]] <- ceiling(chr.lens$length.bp / bin.len)
  
  write.table(x=chr.lens, file=chrlen.filepath, 
              sep="\t", quote=F, append=F, row.names=F, col.names=T)
  
}

# rm(list=ls()); gc()
