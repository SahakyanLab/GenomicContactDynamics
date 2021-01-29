################################################################################
# Filter COSMIC non-coding mutation data
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
CosmicNCVPath = paste0(data.dir, "/cosmic/CosmicNCV.tsv")
out.dir = paste0(wk.dir, "/out_filter")
### OTHER SETTINGS #############################################################
LO.chain = "hg38ToHg19"
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
ncv.df <- fread(file=CosmicNCVPath, header=TRUE) 
if( any(ncv.df$GRCh!=38) ){ stop("Not all are hg38.") }

#-------------------Transform coordinates to 1-based and liftover to hg19
# Extract chr, start, end coordinate from genome position
tmp <- strsplit(x=ncv.df$`genome position`, split=":|-")
rm(ncv.df); gc()
tmp <- lapply(X=tmp, FUN=as.numeric)
tmp <- as.data.frame(do.call("rbind", tmp))
dimnames(tmp)[[2]] <- c("chr", "start", "end")
tmp$chr[tmp$chr==24] <- "Y"
tmp$chr[tmp$chr=="23"] <- "X"
if( any(!tmp$chr%in%c(1:22, "X", "Y")) ){ print("Weird chr.", quote=FALSE) }
print(unique(tmp$chr), quote=FALSE)
tmp$chr <- paste0("chr", tmp$chr)

# Convert 0-based coordinates to 1-based
len <- tmp$end-tmp$start
if( any(len<0) ){ print("Negative length.") }
tmp$start[len>0] <- tmp$start[len>0]+1
# Liftover, hg38 --> hg19
LO.bed <- liftOveR(conversion=LO.chain,
                   space=tmp$chr, 
                   start=as.numeric( tmp$start ),
                   end=as.numeric( tmp$end ),
                   strand="*",
                   getchain=TRUE,
                   rmchain=FALSE,
                   returnGRangesList=FALSE)
orig.len <- nrow(tmp)
if( any(!(1:orig.len)%in%LO.bed$group) ){ print("1. Not all entries lifted over.") }
rm(tmp)

dup.TF <- duplicated(LO.bed$group)
dup.ind <- LO.bed$group[dup.TF] # Index in tmp (orig)
if( any(dup.ind%in%which(len==0)) ){ print("Some single base lengthen after liftover.") }

# Collapse entries split after liftover (Duplicated groups)
tmp1 <- LO.bed[LO.bed$group%in%dup.ind,]
start <- by(data=tmp1$start, INDICES=tmp1$group, FUN=function(x) paste(x, collapse=";"))
end <- by(data=tmp1$end, INDICES=tmp1$group, FUN=function(x) paste(x, collapse=";"))
rm(tmp1)

LO.bed <- LO.bed[!dup.TF,]; rm(dup.TF)
rownames(LO.bed) <- LO.bed$group
LO.bed[names(start),"start"] <- start
LO.bed[names(end),"end"] <- end

LO.bed <- LO.bed[order(LO.bed$group),]
if( !identical(1:orig.len, LO.bed$group) ){ print("2. Not all entries lifted over.") }

rm(start, end); gc()
#-------------------Reload data and append LO coordinates
ncv.df <- fread(file=CosmicNCVPath, header=TRUE) 
ncv.df <- cbind.data.frame(index=1:nrow(ncv.df), 
                           LO.bed[, c("seqnames", "start", "end")], ncv.df)
rm(LO.bed); gc()
ncv.df$"GRCh" <- 37
ncv.df <- ncv.df[-(unique(dup.ind)),]
ncv.df$start <- as.numeric(ncv.df$start)
ncv.df$end <- as.numeric(ncv.df$end)
#-------------------Filter
incl.TF <- ncv.df$`Mutation somatic status`=="Confirmed somatic variant" &
           nchar(ncv.df$WT_SEQ)==1 & nchar(ncv.df$MUT_SEQ)==1 &
           ncv.df$SNP=="n" & #ncv.df$FATHMM_MKL_NON_CODING_SCORE>=0.7 &
           (ncv.df$end-ncv.df$start)==0 & ncv.df$seqnames%in%paste0("chr", c(1:22))
ncv.df <- ncv.df[incl.TF,] 
print(nrow(ncv.df), " final mutations.", quote=FALSE)
rm(incl.TF); gc()
print(paste0(nrow(ncv.df), " mutations after filtering."), quote=FALSE)
save(ncv.df, file=paste0(out.dir, "/CosmicNCV_hg37_final.RData"))

# rm(list=ls()); gc()
