################################################################################
# Function to get complementarity score in terms of kmer or alignment (edit distance)
# Parallel execution by chromosome
# DEPENDENCIES:
#-library(foreach)
#-library(doParallel)
#-library(itertools)
#-library(compiler)
#-library(Biostrings)
#-source(paste0(lib.TrantoRextr, "/UTIL_readLinesFast.R"))
#-source(paste0(lib.TrantoRextr, "/GEN_readfasta.R"))
#-source(paste0(lib.TrantoRextr, "/GEN_loadGenome.R")) 
#-source(paste0(lib.TrantoRextr, "/GEN_getGenomicSeq.R"))
#-source(paste0(wk.dir, "/lib/getComplementarityScore.R"))

#-if(type=="align"){
#-  library(Rcpp)
#-  sourceCpp(file=paste0(wk.dir, "/lib/edlibNW.cpp"))
#-  source(paste0(wk.dir, "/lib/getComplementarityAlign.R"))
#-} else if(type=="kmer"){
#-  source(paste0(wk.dir, "/lib/getComplementarityKmer.R"))
#-} else {
#-  stop("Invalid input for type.")
#-}
################################################################################
################################################################################
getComplementarity <- function(
  # both
  lib.TrantoRextr = paste0(lib, "/TrantoRextr"),
  out.dir = paste0(wk.dir, "/out_constraints"),
  persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc"),
  # File with chromosome lengths (use right genome build), Columns: chromosome-length.bp
  chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt"),
  # align
  genome.dir = paste0(data.dir, "/human_genome_unmasked_37.73"),
  # kmer
  gfreepar.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc"),
  binkmer.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/binkmer_divLen_all"),
  # both
  gcb = "min2Mb",
  chr.v = paste0("chr", c(22:1, "X")),
  bin.len = 4e4L,
  kmer.len = 7L,
  type = "kmer", # "kmer" | "align"
  # For type=align, nCPU based on number of chunks
  # For type=kmer, nCPU based on number of contacts, ~30G for chr1
  # chr21 - align - 368511 good contacts - 30G - 2 days
  nCPU = 3L, # chr1 - 4L, chr21 - 2L
  allij = TRUE,
  # If allij = FALSE, use PERSIST.MX, filter by cell/tissuet type?
  ct = NULL, 
  affix.persist = "",
  affix.binkmer = "", # paste0(chr, "_Hyb", kmer.len, "_", gcb, affix)
  affix.out = "",
  genome.prefix="Homo_sapiens.GRCh37.73.dna.chromosome.",
  fastafile.ending=".fa",
  # align
  numChunks = 30L, # chr1 - 32L, chr21 - 2L
  # kmer
  gfreeparfile = paste0(gfreepar.dir, "/Gfree_", kmer.len, "mer.par")
){
  
  # Chromosome length file
  chrLen.df <- read.table(file=chrLenfile, as.is=TRUE, header=TRUE)
  
  for(chr in chr.v){
    
    print(paste0(chr, "..."), quote=FALSE)
    
    affix.persist <- ifelse(is.null(affix.persist), "", affix.persist)
    affix.binkmer <- ifelse(is.null(affix.binkmer), "", affix.binkmer)
    affix.out <- ifelse(is.null(affix.out), "", affix.out)
    
    start.time <- Sys.time()
    
    # Load from PERSIST.MX to get Cp data; rowname of contact in PERSIST.MX$hits 
    # corresponds to the index of that contact in a matrix of all ij contacts for a
    # given chr obtained using expand.grid() with self and duplicated pairs not yet 
    # removed. The rowname therefore can be used to assign the Cp value of some 
    # contacts in the all ij contact matrix.
    #load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, affix.persist, ".RData"))
    
    # Make mx of contacts
    if(allij){
      
      print("All ij contacts, where i<j...", quote=FALSE)
      
      tot.bin <- ceiling(chrLen.df[chrLen.df$chromosome==chr, "length.bp"]/bin.len)
      ubins <- 1:tot.bin
      
      # Create matrix of all pairs 
      contact.mx <- data.matrix(expand.grid(ubins, ubins))
      dimnames(contact.mx)[[1]] <- 1:nrow(contact.mx)
      rm(ubins, tot.bin)
      
      # Add Cp using the rownames of PERSIST.MX$hits
      contact.mx <- cbind(contact.mx, Cp=NA)
      #contact.mx[rownames(PERSIST.MX$hits), "Cp"] <- PERSIST.MX$ntis
      
      # To remove self pairs and duplicates like {1,2}(keep) and {2,1}
      contact.mx <- contact.mx[contact.mx[,2]-contact.mx[,1]>0,]
      
    } else {
      
      print("ij contacts from PERSIST.MX...", quote=FALSE)
      
      if( ct%in%c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
                  "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC") ){
        
        print("Cell/tissue filtering.", quote=FALSE)
    
        contact.mx <- cbind( data.matrix( PERSIST.MX$hits[, c("i", "j")] ),
                             Cp=PERSIST.MX$ntis )
        contact.mx[PERSIST.MX$hits[,ct]==0,] <- c(NA, NA, NA)
        
      } else if(ct=="hg19"){
        
        print("No Cell/tissue filtering.", quote=FALSE)
        
        contact.mx <- cbind( data.matrix( PERSIST.MX$hits[, c("i", "j")] ),
                             Cp=PERSIST.MX$ntis )
      } else {
        stop("Invalid ct.")
      }
      
    }
    
    #rm(PERSIST.MX); gc()
    
    if(type=="kmer"){
      # Load BINKMER.MX
      load(file=paste0(binkmer.dir, "/", chr, "_BinKmer", kmer.len, affix.binkmer, ".RData"))
      # Normalisa BINKMER.MX to bin length
      BINKMER.MX[,5:ncol(BINKMER.MX)] <- BINKMER.MX[,-(1:4)]/BINKMER.MX[,4] 
    }
    
    # Calculate complementarity score
    CII.MX <- cbind(contact.mx,
                    getComplementarityScore(
                      chr=chr,
                      contact.mx=contact.mx[,1:2],
                      binlength=bin.len,
                      type=type,
                      # For type=align, nCPU based on number of chunks
                      # For type=kmer, nCPU based on number of contacts 
                      nCPU=nCPU,
                      
                      # Alignment-specific arguments
                      PATH.genome=genome.dir,
                      genome.prefix=genome.prefix, 
                      fastafile.ending=fastafile.ending,
                      silent=TRUE,
                      # Number of chunks to split the contacts for alignment
                      numChunks=numChunks,
                      
                      # K-mer-method-specific arguments
                      out.dir=out.dir,
                      out.id=paste0(chr, "_Hyb", kmer.len, "_", gcb, affix.out),
                      binkmer.mx=BINKMER.MX,
                      gfreeparfile=gfreeparfile,
                      # Kmer length
                      kmerlength=kmer.len,
                      saveOut=TRUE
                    )
    )
    rm(contact.mx); gc()
    
    dimnames(CII.MX)[[2]] <- c("i", "j", "Cp", "C||")
    save(CII.MX, file=paste0(out.dir, "/", chr, "_", type, "_", gcb, affix.out, ".RData"))
    
    rm(CII.MX); gc()
    
    end.time <- Sys.time()
    print(end.time-start.time, quote=FALSE)
    
    print(paste0(chr, " done!"), quote=FALSE)
    
  }
  
}
################################################################################
suppressWarnings(suppressPackageStartupMessages(library(compiler)))
getComplementarity <- cmpfun(getComplementarity, options=list(suppressUndefined=TRUE))
################################################################################
