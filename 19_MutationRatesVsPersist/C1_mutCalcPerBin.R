################################################################################
# Take mutation file, assign mutated sites to bin and then calculate per bin:
# a. Tmut - number of mutations 
# b. Nmsite - number of mutated sites
# c. TmutDIVNmsite - Tmut/Nmsite
# d. Nmsitenorm - Nmsite/number of wildtype base in the bin 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
src.dir = paste0(wk.dir, "/out_filter")
basecont.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binBaseContent")
out.dir = paste0(wk.dir, "/out_mutCalcPerBin")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
src.id = "hg38ToHg19" # "Hg19" | "hg38ToHg19"
bin.len = 40000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrlen.df <- read.delim(file=chrLenfile, header=TRUE)

#load(file=paste0(src.dir, "/CosmicNCV_", src.id, "_final_2000.RData"))
load(file=paste0(src.dir, "/CosmicNCV_", src.id, "_final.RData"))
ncv.df <- ncv.df[,c("ID_SAMPLE", "chr", "start", "MUT")]
#-------------------Assign mutated site to bin
ncv.df$bin <- ceiling(ncv.df$start/bin.len) 

chr.v <- unique(ncv.df$chr)

mut.v <- c("All", unique(ncv.df$MUT))

for(mut in mut.v){
  
  if(mut=="All"){
    mut.TF <- rep(TRUE, times=nrow(ncv.df))
    WT_SEQ <- "numUMChar"
  } else if( mut%in%c("A>C", "A>G", "A>T", "C>A", "C>G", "C>T") ){
    mut.TF <- ncv.df$MUT==mut
    WT_SEQ <- strsplit(x=mut, split=">", fixed=TRUE)[[1]][1]
  } else {
    stop("Invalid mutation type notation.")
  }
  
  MUTBIN.DF <- list()
  for(chr in chr.v){
    
    incl.ind <- which(ncv.df$chr==chr & mut.TF)
    if( length(incl.ind)==0 ){
      next
      print(paste0(chr, " ", mut, ":Skipped."))
    }
      
    chr.bin <- ceiling((chrlen.df$length.bp[chrlen.df$chromosome==chr])/bin.len)
    if( any(ncv.df[incl.ind,"bin"]>chr.bin) ){
      stop(paste0(chr, " ", mut, ":Bin exceeding chr."))
    }
    
    # Total mutation per bin
    Tmut <- table(ncv.df[incl.ind,"bin"])
    
    # Total mutation per bin divided by number of mutated sites
    x <- unique(paste0(ncv.df$start[incl.ind], "..", ncv.df$bin[incl.ind]))
    x <- unlist(lapply(X=strsplit(x=x, split="..", fixed=TRUE), FUN=function(x)x[2]))
    Nmsite <- table(x) # Number of mutated sites
    bin.srt <- as.character(sort(as.numeric(names(Tmut)), decreasing=FALSE))
    if( !all(bin.srt%in%names(Nmsite)) ){
      stop(paste0(chr, " ", mut, ": Checkpoint 1."))
    }
    TmutDIVNmsite <- Tmut[bin.srt]/Nmsite[bin.srt]
    
    # Number of mutated sites normalised to number of original/WT base per bin
    load(file=paste0(basecont.dir, "/", chr, "_BinKmer1.RData"))
    dimnames(BINKMER.MX)[[1]] <- BINKMER.MX[,"bins"]
    numWTSEQ <- BINKMER.MX[bin.srt,WT_SEQ]
    Nmsitenorm <- Nmsite[bin.srt]/numWTSEQ
    rm(BINKMER.MX, numWTSEQ); gc()
    
    MUTBIN.DF[[chr]] <- cbind.data.frame(chr=chr, bin=as.numeric(bin.srt), 
                                         Tmut=as.numeric(Tmut[bin.srt]), 
                                         Nmsite=as.numeric(Nmsite[bin.srt]),
                                         TmutDIVNmsite=as.numeric(TmutDIVNmsite[bin.srt]),
                                         Nmsitenorm=as.numeric(Nmsitenorm),
                                         stringsAsFactors=FALSE)
    print(paste0(chr, " done!"), quote=FALSE)
    
  } # chr.v for loop end
  
  MUTBIN.DF <- do.call("rbind", MUTBIN.DF)
  rownames(MUTBIN.DF) <- NULL
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=TRUE)
  save(MUTBIN.DF, file=paste0(out.dir, "/", mut.id, "_", src.id, "_mutCalcPerBin.RData"))
  
}

# rm(list=ls()); gc()

