################################################################################
# Liftover and apply series of filtering to get final cosmic non-variant dataset.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
genome = 37 # 37 | 38
CosmicNCVPath = paste0(data.dir, "/cosmic/GRCh", genome, "/CosmicNCV.GRCh", 
                       genome, ".tsv")
# 1-based, hg19, beds of considered coding areas
bed.dir  = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/regions_anno/out_extractRegions")
out.dir = paste0(wk.dir, "/out_filter")
### OTHER SETTINGS #############################################################
LO.chain = "Hg19" # "Hg19" | "hg38ToHg19"
sample.len = 100000
filter3Only = FALSE
bed.header = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(rtracklayer)
library(GenomicRanges)
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
source(paste0(lib, "/TrantoR_liftOver/GEN_liftOveR.R"))
source(paste0(lib, "/TrantoR_liftOver/GEN_liftOverLoadChain.R"))
source(paste0(wk.dir, "/lib/applyFilter1.R"))
source(paste0(wk.dir, "/lib/applyFilter2.R"))
source(paste0(wk.dir, "/lib/applyFilter3.R"))
source(paste0(lib, "/whichOverlapToBed.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(filter3Only==FALSE){
  
  ncv.df <- fread(file=CosmicNCVPath, header=TRUE, data.table=FALSE, stringsAsFactors=FALSE) 
  if( any(ncv.df$GRCh!=genome) ){ stop(paste0("Not all are GRCh", genome)) }
  
  #-------------------GET WGS mutations
  ncv.df <- ncv.df[ncv.df$Whole_Genome_Reseq=="y" & ncv.df$Whole_Exome=="n",]
  
  #-------------------Apply applyFilter1()
  
  filter1.TF <- applyFilter1(ncv.df=ncv.df)
  ncv.df <- ncv.df[filter1.TF$SAMP.TF & filter1.TF$NOTREDUNDANT.TF,]
  rm(filter1.TF)
  
  #-------------------Apply applyFilter2()
  
  filter2.TF <- applyFilter2(ncv.df=ncv.df)
  ncv.df <- ncv.df[!(filter2.TF$GMI.drop.TF | filter2.TF$POS.drop.TF),]
  rm(filter2.TF)
 
  gc()
  
  save(ncv.df, file=paste0(out.dir, "/CosmicNCV_", LO.chain, "_WGSfiltered1And2.RData"))
  
  #-------------------Extract chr, start, end coordinate from genome position
  tmp <- strsplit(x=ncv.df$`genome position`, split=":|-")
  rm(ncv.df); gc()
  tmp <- lapply(X=tmp, FUN=function(x){
    
    if(length(x)!=3){
      stop("Problematic genome position")
    }
    return(as.numeric(x))
    
  })
  tmp <- as.data.frame(do.call("rbind", tmp))
  
  # Fix numeric chr notation
  dimnames(tmp)[[2]] <- c("chr", "start", "end")
  tmp$chr <- as.character(tmp$chr)
  tmp$chr[tmp$chr=="23"] <- "X"
  tmp$chr[tmp$chr=="24"] <- "Y"
  tmp$chr[tmp$chr=="25"] <- "MT" # Not sure about this
  if( any(!tmp$chr%in%c(1:22, "X", "Y", "MT")) ){ stop("Weird chr.", quote=FALSE) }
  tmp$chr <- paste0("chr", tmp$chr)
  
  len <- tmp$end-tmp$start
  if( any(len<0) ){ print("One-based but some have negative length.") }
  orig.len <- nrow(tmp)
  print(paste0(orig.len, " WGS mutations after filters 1 and 2..."), quote=FALSE)
  
  #-------------------If required, liftover to hg19/GRCh37
  if(genome!=37){
    
    print(paste0("Liftover coordinates ", LO.chain, "..."), quote=FALSE)
    
    # Liftover; recognises chrX and chrY (not numeric equivalent)
    tmp <- liftOveR(conversion=LO.chain,
                    space=tmp$chr, 
                    start=as.numeric(tmp$start),
                    end=as.numeric(tmp$end),
                    strand="*",
                    getchain=TRUE,
                    rmchain=FALSE,
                    returnGRangesList=FALSE)
    data.table::setnames(x=tmp, old="seqnames", new="chr")
    tmp$chr <- as.character(tmp$chr)
    tmp$strand <- as.character(tmp$strand)
    if( is.factor(tmp$group) ){ "Group column is a factor." }
    tmp$group <- as.numeric(as.character(tmp$group))
    if( any(!(1:orig.len)%in%tmp$group) ){
      print("1. Not all entries lifted over.")
    }
    
    # Index of mutations in tmp (orig) that are duplicated in tmp (LO), meaning split after liftover
    dupgroup.v <- tmp$group[duplicated(tmp$group)]
    if( any(dupgroup.v%in%which(len==0)) ){
      print("Some single base lengthen after liftover.")
    }
    # Remove mutations that were split after liftover
    tmp <- tmp[!as.numeric(tmp$group)%in%as.numeric(dupgroup.v),]
    rm(dupgroup.v)
    if( any(duplicated(tmp$group)) ){ stop("Duplicated mutation.") }
    
  } else {
    
    print("GRCh37 coordinates. No need to liftOver...", quote=FALSE)
    LO.chain <- "Hg19"
    tmp$group <- 1:orig.len
    
  }
  
  #-------------------Reload data and add LO coordinates
  load(paste0(out.dir, "/CosmicNCV_", LO.chain, "_WGSfiltered1And2.RData"))
  ncv.df$start <- ncv.df$end <- NA
  ncv.df$GRCh <- 37
  tmp$group <- as.numeric(tmp$group)
  ncv.df$chr[tmp$group] <- tmp$chr
  ncv.df$start[tmp$group] <- as.numeric(tmp$start)
  ncv.df$end[tmp$group] <- as.numeric(tmp$end)
  rm(tmp); gc()
  save(ncv.df, file=paste0(out.dir, "/CosmicNCV_", LO.chain, "_WGSfiltered1And2.RData"))
  
} else {
  load(paste0(out.dir, "/CosmicNCV_", LO.chain, "_WGSfiltered1And2.RData"))
}

#-------------------Rest

#Apply applyFilter3()
ncv.df <- ncv.df[applyFilter3(ncv.df=ncv.df)$incl.TF,]

ncv.df$MUT <- paste0(ncv.df$WT_SEQ, ">", ncv.df$MUT_SEQ)
# COSMIC convention
ncv.df$MUT[ncv.df$MUT=="A>C"] <- "T>G"
ncv.df$MUT[ncv.df$MUT=="A>G"] <- "T>C"
ncv.df$MUT[ncv.df$MUT=="A>T"] <- "T>A"
ncv.df$MUT[ncv.df$MUT=="G>T"] <- "C>A"
ncv.df$MUT[ncv.df$MUT=="G>C"] <- "C>G"
ncv.df$MUT[ncv.df$MUT=="G>A"] <- "C>T"

if( any(!ncv.df$MUT%in%c("T>G","T>C","T>A","C>A","C>G","C>T")) ){
  stop("Problematic WT_SEQ/MUT_SEQ.")
}

# Classify mutations into coding (C) and non-coding (NC)
ncv.df$Coding <- whichOverlapToBed(bed.dir=bed.dir, header=bed.header, seqnames=ncv.df$chr, 
                                   start=ncv.df$start, end=ncv.df$end, 
                                   strand=rep("*", times=length(ncv.df$start))
                                   )

# Final dataset after series of filtering
save(ncv.df, file=paste0(out.dir, "/CosmicNCV_", LO.chain, "_final.RData"))
print(paste0(nrow(ncv.df), " mutations after filtering."), quote=FALSE)

# Extract sample
set.seed(123)
ncv.df <- ncv.df[sample(x=1:nrow(ncv.df), size=sample.len),]
sample.len <- format(x=nrow(ncv.df), scientific=FALSE)
save(ncv.df, file=paste0(out.dir, "/CosmicNCV_", LO.chain, "_final_", sample.len, ".RData"))

# rm(list=ls()); gc()
