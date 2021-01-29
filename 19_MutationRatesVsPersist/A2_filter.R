################################################################################
# Apply filtering based on checks done in the previous script, liftover to hg19 
# and conversion to 1-based coordinate system to get the final non-coding
# mutation working dataset. 
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
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
genome = 38 # 37 | 38
CosmicNCVPath = paste0(data.dir, "/cosmic/GRCh", genome, "/CosmicNCV.GRCh", 
                       genome, ".tsv")
out.dir = paste0(wk.dir, "/out_filter")
### OTHER SETTINGS #############################################################
LO.chain = "hg38ToHg19" # "Hg19" | "hg38ToHg19"
filterOnly = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(rtracklayer)
library(GenomicRanges)
source(paste0(lib, "/TrantoR_liftOver/GEN_liftOveR.R"))
source(paste0(lib, "/TrantoR_liftOver/GEN_liftOverLoadChain.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(filterOnly==FALSE){
  
  ncv.df <- fread(file=CosmicNCVPath, header=TRUE, data.table=TRUE) 
  if( any(ncv.df$GRCh!=genome) ){ stop(paste0("Not all are GRCh", genome)) }
  
  # Check for NAs per column
  tmp <- list()
  for( i in 1:ncol(ncv.df) ){
    colnme <- colnames(ncv.df)[i]
    if( any(is.na(ncv.df[[colnme]])) ){
      tmp[[colnme]] <- colnme
    }
    rm(colnme)
  }
  print("Columns with NAs...", quote=FALSE)
  print(names(tmp), quote=FALSE)
  
  #-------------------Extract chr, start, end coordinate from genome position
  tmp <- strsplit(x=ncv.df$`genome position`, split=":|-")
  rm(ncv.df); gc()
  tmp <- lapply(X=tmp, FUN=as.numeric)
  tmp <- as.data.frame(do.call("rbind", tmp))
  
  # Fix numeric chr notation
  dimnames(tmp)[[2]] <- c("chr", "start", "end")
  tmp$chr[tmp$chr==23] <- "X"
  tmp$chr[tmp$chr=="24"] <- "Y"
  tmp$chr[tmp$chr=="25"] <- "MT" # Not sure about this
  if( any(!tmp$chr%in%c(1:22, "X", "Y", "MT")) ){ print("Weird chr.", quote=FALSE) }
  print(unique(sort(tmp$chr)), quote=FALSE)
  tmp$chr <- paste0("chr", tmp$chr)
  
  # Convert 0-based coordinates to 1-based
  len <- tmp$end-tmp$start
  if( any(len<0) ){ print("Negative length.") }
  tmp$start[len>0] <- tmp$start[len>0]+1
  orig.len <- nrow(tmp)
  print(paste0(orig.len, " original mutations"), quote=FALSE)
  
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
    if( any(!(1:orig.len)%in%as.numeric(tmp$group)) ){ print("1. Not all entries lifted over.") }
    
    dup.TF <- duplicated(tmp$group)
    # Index of mutations in tmp (orig) that are duplicated in tmp (LO), meaning split after liftover
    dup.ind <- tmp$group[dup.TF]
    if( any(dup.ind%in%which(len==0)) ){ print("Some single base lengthen after liftover.") }
    # Remove mutations that were split after liftover
    tmp <- tmp[!tmp$group%in%dup.ind,]
    rm(dup.ind, dup.TF)
    if( any(duplicated(tmp$group)) ){ stop("Duplicated mutation.") }
    
  } else {
    
    print("GRCh37 coordinates. No need to liftOver...", quote=FALSE)
    LO.chain <- "Hg19"
    tmp$group <- 1:orig.len
    
  }
  
  #-------------------Reload data and add LO coordinates
  ncv.df <- fread(file=CosmicNCVPath, header=TRUE, data.table=FALSE) 
  ncv.df$chr <- ncv.df$start <- ncv.df$end <- NA
  tmp$group <- as.numeric(tmp$group)
  ncv.df$chr[tmp$group] <- tmp$chr
  ncv.df$start[tmp$group] <- as.numeric(tmp$start)
  ncv.df$end[tmp$group] <- as.numeric(tmp$end)
  rm(tmp); gc()
  save(ncv.df, file=paste0(out.dir, "/CosmicNCV_", LO.chain, ".RData"))
  print(paste0(nrow(ncv.df), "/", orig.len, " mutations."), quote=FALSE)
  
} else {
  load(file=paste0(out.dir, "/CosmicNCV_", LO.chain, ".RData"))
}

#-------------------Filter

# Remove entries at same position and potentially from same sample; taking out
# this filter only adds ~100K entries
#ncv.df$custom.id <- paste(ncv.df$`Primary site`, ncv.df$`Site subtype 1`, 
#                          ncv.df$`Site subtype 2`, ncv.df$`Site subtype 3`,
#                          ncv.df$`Primary histology`, ncv.df$`Histology subtype 1`,
#                          ncv.df$`Histology subtype 2`, ncv.df$`Histology subtype 3`)

for( x in c("ID_SAMPLE", "Sample name") ){ #, "custom.id") ){
  
  drp <- paste0(ncv.df[[x]], ncv.df$`genome position`)
  drp <- duplicated(drp) | duplicated(drp, fromLast=TRUE)
  if(x=="ID_SAMPLE"){
    SiDdrop.TF <- drp
  } else {
    SiDdrop.TF <- SiDdrop.TF | drp
  }
  rm(drp)
  
}

# Remove entries with same GENOMIC_MUTATION_ID but different position
dupGid <- ncv.df$GENOMIC_MUTATION_ID[duplicated(ncv.df$GENOMIC_MUTATION_ID)]
dupGid.ind <- which(ncv.df$GENOMIC_MUTATION_ID%in%dupGid)
same.TF <- by(data=ncv.df$`genome position`[dupGid.ind],
              INDICES=ncv.df$GENOMIC_MUTATION_ID[dupGid.ind],
              FUN=function(x) length(unique(as.character(x)))==1
)
GidIncl.TF <- rep(TRUE, times=length(SiDdrop.TF))
GidIncl.TF[ dupGid.ind[!same.TF] ] <- FALSE

# Rest of filters
base.v <- c("A", "T", "C", "G")
base.v <- c(base.v, tolower(base.v))
incl.TF <- ncv.df$`Mutation somatic status`=="Confirmed somatic variant" &
  # Removing SNP filter only adds ~60K entries
  #ncv.df$SNP=="n" & #ncv.df$FATHMM_MKL_NON_CODING_SCORE>=0.7 & most have missing scores
  ncv.df$WT_SEQ%in%base.v & ncv.df$MUT_SEQ%in%base.v & !SiDdrop.TF & GidIncl.TF &
  (ncv.df$end-ncv.df$start)==0 & ncv.df$chr%in%paste0("chr", c(1:22)) &
  ncv.df$Whole_Genome_Reseq=="y" & ncv.df$Whole_Exome=="n"
ncv.df <- ncv.df[incl.TF,] 
rm(incl.TF); gc()

ncv.df$MUT <- paste0(ncv.df$WT_SEQ, ">", ncv.df$MUT_SEQ)
ncv.df$MUT[ncv.df$MUT=="T>G"] <- "A>C"
ncv.df$MUT[ncv.df$MUT=="T>C"] <- "A>G"
ncv.df$MUT[ncv.df$MUT=="T>A"] <- "A>T"
ncv.df$MUT[ncv.df$MUT=="G>T"] <- "C>A"
ncv.df$MUT[ncv.df$MUT=="G>C"] <- "C>G"
ncv.df$MUT[ncv.df$MUT=="G>A"] <- "C>T"

if( any(!ncv.df$MUT%in%c("A>C","A>G","A>T","C>A","C>G","C>T")) ){
  stop("Problematic WT_SEQ/MUT_SEQ.")
}

# Dataset after filtering
save(ncv.df, file=paste0(out.dir, "/CosmicNCV_", LO.chain, "_final.RData"))
print(paste0(nrow(ncv.df), " mutations after filtering."), quote=FALSE)

# Extract sample
ncv.df <- ncv.df[1:2000,]
save(ncv.df, file=paste0(out.dir, "/CosmicNCV_", LO.chain, "_final_2000.RData"))

# rm(list=ls()); gc()
