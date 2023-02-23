################################################################################
# Using all mutCalcPerBin.RData, generate per chr and per mut metric/calc an ij.mut 
# object containing consensus value of mut calc per contact based on ij.fnx(na.rm=F)
# (order matches with PERSIST.MX). 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/19_MutationRatesVsPersist")
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
src.dir = paste0(wk.dir, "/out_mutCalcPerBin")
out.dir = paste0(wk.dir, "/out_mut_contact_Cp_plotdata")
### OTHER SETTINGS #############################################################
chr = "chrCHRREPLACE" #chrs = paste0("chr", c(21:22))
nCPU = 5 # src.files
persist.id = "Persist_min2Mb"
src.file.pattern = c("donor_centric_PCAWG_Hg19", "sigEperclimits_1000rawInf_" , 
                     "mutCalcPerBin.RData") # case-sensitive
ij.fnx = "mean" # "median" # IJFNXREPLACE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Get list of files satisfying all patterns

all.files <- list.files(path=src.dir)
is.src <- sapply(src.file.pattern, simplify=T, FUN=function(patt){
  grepl(patt, all.files, ignore.case=F, fixed=T)
})
src.files <- all.files[ rowSums(is.src) == length(src.file.pattern) ]

# Load contact data

load(paste0(persist.dir, "/", chr, "_", persist.id, ".RData"))
ij.mx <- cbind(i=PERSIST.MX$hits$i, j=PERSIST.MX$hits$j, Cp=PERSIST.MX$ntis)
rm(PERSIST.MX)

# Generate IJ.MUT per src.file 

src.files.len <- length(src.files)
#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:src.files.len, chunks=nCPU), .inorder=F
               
) %op% {
  
  for(i in itr){
    
    src.file <- src.files[[i]]
    
    # Load mut metric bin data
    
    load(paste0(src.dir, "/", src.file))
    
    mut.calcs <- setdiff( colnames(MUTBIN.DF), c("chr", "bin") )
    
    MUTBIN.DF <- MUTBIN.DF[MUTBIN.DF$chr %in% chr,]
    
    # Note if there's no data for chr
    
    if( length(MUTBIN.DF$chr) == 0 ){
      message(paste0(src.file, " ", chr, ": No data, skipped."))
      rm(MUTBIN.DF)
      next
    }
    
    out.id <- paste0("ijfnx", ij.fnx, "_", gsub("mutCalcPerBin.", "ijmut.", src.file, fixed=T))
    
    #toExport <- c("chrs", "persist.dir", "persist.id", "MUTBIN.DF", "mut.calcs", "out.id",
    #              "ij.fnx", "out.dir")
    # #### PARALLEL EXECUTION #########
    # foreach(itr=isplitVector(1:chrs.len, chunks=nCPU), .inorder=F,
    #         .export=toExport, .noexport=ls()[!ls()%in%toExport]
    #               
    # ) %op% {
    #   
    #   for(i in itr){
    #     
    #     chr <- chrs[[i]]
    
    # From MUTBIN.DF, make matrix of bin mut.calcs based on order of i and j bins in contact matrix
    
    is.chr <- MUTBIN.DF$chr == chr
    
    ind.MUTBINDF <- match(ij.mx[,"i"], table=MUTBIN.DF$bin[is.chr])
    i.mut.mx <- data.matrix( MUTBIN.DF[is.chr, mut.calcs][ind.MUTBINDF,] )
    
    ind.MUTBINDF <- match(ij.mx[,"j"], table=MUTBIN.DF$bin[is.chr])
    j.mut.mx <- data.matrix( MUTBIN.DF[is.chr, mut.calcs][ind.MUTBINDF,] )
    
    if( !identical(dim(i.mut.mx), dim(j.mut.mx)) ){
      rm(MUTBIN.DF)
      stop("i.mut.mx and j.mut.mx differ in dimension.")
    }
    
    rownames(i.mut.mx) <- rownames(j.mut.mx) <- NULL
    rm(ind.MUTBINDF)
    
    # Apply ij.funx from generated i and j matrices to derive consensus metric per contact
    
    if(ij.fnx == "mean"){
      ij.mut.mx <- (i.mut.mx + j.mut.mx) / 2
    } else {
      
      ij.mut.array <- simplify2array(list(i.mut.mx, j.mut.mx))
      eval(parse(text=paste0(
        'ij.mut.mx <- apply(ij.mut.array, MARGIN=c(1,2), FUN=', ij.fnx, ', na.rm=F)'
      )))
      
    }
    
    rm(i.mut.mx, j.mut.mx)
    
    # Add Cp to ij.mut.mx
    ij.mut.mx <- cbind(ij.mut.mx, Cp=ij.mx[,"Cp"])
    #rm(ij.mx)
    
    # Save IJ.MUT per mut metric
    
    for(calc in mut.calcs){
      
      IJ.MUT <- cbind(value=ij.mut.mx[,calc], Cp=ij.mut.mx[,"Cp"])
      save(IJ.MUT, file=paste0(out.dir, "/", chr, "_", calc, "_", out.id))
      rm(IJ.MUT)
      
    } # mut.calcs for loop end
    
    rm(ij.mut.mx)
    
    #  }
    
    # } # chrs.len foreach loop end
    # ### END OF PARALLEL EXECUTION ###
    
    rm(MUTBIN.DF)
    
    message(paste0(src.file, " done!"))
    
  } # src.files for loop end
  
}
### END OF PARALLEL EXECUTION ###
  
  
# rm(list=ls()); gc()