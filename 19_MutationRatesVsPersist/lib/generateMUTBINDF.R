################################################################################
# Generate MUTBIN.DF (chr combined) per combination of mutation type,
# signature/s of interest, signature exposure threshold and mutation location.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(compiler)
# source(paste0(wk.dir, "/lib/selectSamplesBasedOnSigExposure.R"))
# source(paste0(wk.dir, "/lib/makeMUTBINDFperChrPerMUT.R"))
### FUNCTION ###################################################################
generateMUTBINDF <- function(ncv.df, data.id, src.id, out.dir, basecont.dir, 
                             basecont.affix, bin.len, chrlen.df, mut, sigExposure.df, 
                             SIG, sigEpLim, loc
){
  
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=T)
  
  #-------------------FILTERS:
  
  # Samples/donors based on signature exposure to signature/s of interest
  sigEpLim.id <- sigEpLim
  sigEpLim <- strsplit(x=sigEpLim, split="_", fixed=T)[[1]]
  
  SIG.id <- SIG
  SIG <- strsplit(x=SIG, split="_", fixed=T)[[1]]
  
  if( all(SIG%in%colnames(sigExposure.df)) & sigEpLim.id!="nosampfilter" ){
    
    sample.v <- selectSamplesBasedOnSigExposure(sigExposure.df=sigExposure.df, sig.v=SIG, 
                                                sigEperc.limits=sort(as.numeric(sigEpLim), 
                                                                     decreasing=F)
    )
    
  } else if(sigEpLim.id=="nosampfilter"){
    print('sigEpLim.id="nosampfilter"; Using all samples...', quote=F)
  } else {
    stop("Invalid input signatures...", quote=F)
  }
  
  if( sigEpLim.id=="nosampfilter" & !exists(x="sample.v") ){
    sample.v <- unique(ncv.df$ID_SAMPLE)
  }
  
  # Location  
  loc.id <- loc
  loc <- strsplit(x=loc, split="_", fixed=T)[[1]]
  
  # Mutation type
  if(mut=="All"){
    incl.TF <- ncv.df$location%in%loc & ncv.df$ID_SAMPLE%in%sample.v
  } else if( mut%in% c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G") ){
    incl.TF <- ncv.df$location%in%loc & ncv.df$ID_SAMPLE%in%sample.v & ncv.df$MUT==mut
  } else {
    stop(mut, ": Invalid mutation type notation...")
  }
  
  # Subset mutation data
  ncv.df <- ncv.df[incl.TF,]
  
  if( sum(incl.TF)==0 ){
    
    print( paste0("No mutations left for arguments, ", mut.id, ",", SIG.id, ",",
                  loc, ",", sigEpLim.id), quote=F )
    
    tot.mut <- samp.len <- 0
    Nsamp.mean <- Nsamp.med <- samp.str <- NA
   
  } else {
    
    
    # Get metadata 
    tmp <- table(ncv.df$ID_SAMPLE)
    samp.str <- paste(names(tmp),collapse=";")
    samp.len <- length(tmp)
    Nsamp.mean <- mean(tmp)
    Nsamp.med <- median(tmp)
    
    rm(sample.v, sigEpLim, SIG, loc, incl.TF, tmp)
    
    #-------------------MUTBIN.DF
    
    chr.v <- unique(ncv.df$chr)
    MUTBIN.DF <- list()
    for(chr in chr.v){
      
      load(file=paste0(basecont.dir, "/", chr, "_BinKmer1", basecont.affix, ".RData"))
      chr.len <- chrlen.df$length.bp[chrlen.df$chromosome==chr]
      
      MUTBIN.DF[[chr]] <- makeMUTBINDFperChrPerMUT(ncv.df=ncv.df[ncv.df$chr==chr,], 
                                                   BINKMER.MX=BINKMER.MX, chr=chr, 
                                                   mut=mut, chr.len=chr.len,
                                                   bin.len=bin.len)
      #print(chr)
      
      rm(BINKMER.MX, chr.len); gc()
      
    } # chr.v for loop end
    
    MUTBIN.DF <- do.call("rbind.data.frame", c(MUTBIN.DF, stringsAsFactors=F))
    rownames(MUTBIN.DF) <- NULL
    out.id <- paste0(data.id, "_", src.id, "_", mut.id, "_", SIG.id, "_", loc.id, 
                     "_sigEperclimits_", sigEpLim.id)
    tot.mut <- sum(MUTBIN.DF$Tmut)
    save(MUTBIN.DF, file=paste0(out.dir, "/", out.id, "_mutCalcPerBin.RData"))
    
    print(paste0(out.id, ": MUTBIN.DF generated."), quote=F)
    
  }
  
  print(paste0(samp.len, " samples..."), quote=F)
  
  # Return metadata
  meta <- c(mut.id, SIG.id, loc.id, sigEpLim.id, tot.mut, samp.len, 
            Nsamp.mean, Nsamp.med, samp.str)
  names(meta) <- c("mut.id", "SIG.id", "loc.id", "sigEpLim.id", "Nmut", "Nsamp",
                   "meanNsamp", "medianNsamp", "sample")
  return(meta)
  
}
################################################################################
generateMUTBINDF <- cmpfun(generateMUTBINDF, options=list(suppressUndefined=T))
################################################################################

# rm(list=ls()); gc()