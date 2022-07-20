################################################################################
# Function to mask a or any combination of genomic features and then re-calculate 
# kmer-based discordance score of Hi-C contacts
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(IRanges)
# library(Biostrings)
# library(foreach)
# library(doParallel)
# library(itertools)
# library(compiler)

# source(paste0(lib.TrantoRextr, "/UTIL_readLinesFast.R"))
# source(paste0(lib.TrantoRextr, "/TrantoRextr/GEN_readfasta.R"))
# source(paste0(lib.TrantoRextr, "/GEN_loadGenome.R")) 

# source(paste0(wk.dir, "/lib/maskGenome.R"))

# source(paste0(lib.TrantoRextr, "/GEN_getKmers.R"))                 
# source(paste0(lib.TrantoRextr, "/GEN_getKmerHybridisationGs.R")) 
# source(paste0(lib.TrantoRextr, "/GEN_getKmerCountsPerInterval.R"))  

# source(paste0(wk.dir, "/lib/HiCHybridisation.R"))
# source(paste0(wk.dir, "/lib/HiCHybridPlot.R"))
# source(paste0(wk.dir, "/lib/HiCHybridMaskedGenFeatures.R"))
### FUNCTION ###################################################################
getMaskedKmerCount <- function(
  
  # original BINKMER.MX all bins
  binkmer.dir = "/dir",
  out.dir = "/dir",

  gcb = "min2Mb",
  chr = "chr21", # any of c(1:22,"X","Y","MT")
  kmer.len = 7L,
  HiC.res = 4e4L,
  # Identifier of source data, also used for outputs; _<affix>
  affix = "",
  nCPU = 2L,
  
  PATH.genome = "/Volumes/Data/Database/human_genome_unmasked_37.73",
  genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
  fastafile.ending = ".fa",
  split = FALSE, # sequence splitting while loading fasta
  # If loaded genome exists (hence will not be freshly
  # created by the current loadGenome call, <split> should
  # match the loaded genome split pattern.
  
  # Masking of chromosome skipped if maskingChar or maskbed is NULL
  maskbed = "bed-format table",
  maskingChar = "m" 
  
){

  #-------------------------------------------------------------------------------
  # Load chromosome
  #-------------------------------------------------------------------------------
  chr.id <- strsplit(x=chr, split="chr")[[1]][2]

  loadGenome(PATH.genome=PATH.genome, genome.prefix=genome.prefix,
             fastafile.ending=fastafile.ending,
             chr.id=chr.id, 
             silent = FALSE, split=split,
             remove.other.loads=TRUE)
  
  genome.filename <- paste0(genome.prefix, chr.id, fastafile.ending)
  
  #-------------------------------------------------------------------------------
  # Masked version of chromosome
  #-------------------------------------------------------------------------------
  if( !is.null(maskingChar) | is.null(maskbed) ){
    
    # Mask chromosome
    eval(parse(text=paste0(
      genome.filename, 
      "$seq <- maskGenome(chr=chr, seq=", 
      genome.filename, "$seq, maskbed=maskbed, maskingChar=maskingChar, reduce=TRUE, split=split)"
    )
    ))
    
    # Check if masking worked
    eval(parse(text=paste0(
      "if( grepl(x=", genome.filename, "$seq, pattern=maskingChar) ){ print('Sequence has masked nt.') }"
    )
    ))
    
    print("Genome masking done.", quote=FALSE)
    
  } else {
    
    print("No masking of genome.", quote=FALSE)
    
  }
    
  #-------------------------------------------------------------------------------
  # Count kmers in HiC contact bins (generate masked BINKMER.MX)
  #-------------------------------------------------------------------------------
  # Load unmasked BINKMER.MX (all bins) to get positions of bins
  load(paste0(binkmer.dir, "/", chr, "_BinKmer7.RData"))
  
  ubins <- BINKMER.MX[,"bins"]
  bin.end <- BINKMER.MX[,"endpos"]
  bin.start <- BINKMER.MX[,"startpos"]
  
  rm(BINKMER.MX); gc()

  eval(parse(text=paste0(
               "BINKMER.MX <- getKmerCountsPerInterval(chr=chr, sequence=", genome.filename, "$seq, startpos=bin.start, endpos=bin.end, K=kmer.len, PATH.genome=PATH.genome, genome.prefix=genome.prefix, fastafile.ending=fastafile.ending, silent = FALSE, nCPU=nCPU, maskingChar=maskingChar)"
  )
  ))
  
  # Remove chromosome sequence
  #eval(parse(text=
  #             paste0("rm(", genome.filename, ")")
  #))

  BINKMER.MX <- cbind(bins=ubins, startpos=bin.start, endpos=bin.end, BINKMER.MX)
  rm(ubins, bin.end, bin.start); gc()
  
  save(BINKMER.MX, file=paste0(out.dir, "/", chr, "_BinKmer", kmer.len, affix, ".RData"))
  
}
################################################################################
getMaskedKmerCount <- cmpfun(getMaskedKmerCount, options=list(suppressUndefined=TRUE))
################################################################################