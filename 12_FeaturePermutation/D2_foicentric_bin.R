################################################################################
# Determine what Cp most of the features overlap to
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
binmx.dir = paste0(wk.dir, "/out_bindata")
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/reduced_b1b2b3")
foifile = NULL #paste0(wk.dir, "/foifile/foifile_test")  
out.dir = paste0(wk.dir, "/out_foicentric")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
bin.len = 40000
Cp.v = 1:21
Ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", 
         "SB", "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC",
         "FC", "LC")
out.id = "b1b2b3"
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(RColorBrewer)
library(GenomicRanges)
library(ggplot2)
library(regioneR)
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
source(paste0(lib, "/finaliseFOI.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, " done..."), quote=FALSE)
Cp.v <- sort(Cp.v, decreasing=FALSE)
out.name <- paste0(gcb, "_", out.id)

if(plotOnly==FALSE){
  
  #-------------------List of features and corresponding cell/tissue
  foi.v <- unique(finaliseFOI(foi.dir=foi.dir, foifile=foifile))
  foi.v.len <- length(foi.v)
  foilab.v <- gsub(x=foi.v, pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
  ct.v <- unlist(
    lapply(X=strsplit(x=foi.v, split="ct_|\\_"), FUN=function(x)x[2])
  )
  #-------------------Max Cp of chr bins
  if( is.null(chr.v) ){ chr.v <- paste0("chr", c(1:22, "X")) }
  BIN.MX <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
    load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
    BIN.MX <- cbind.data.frame(chr=chr, BIN.MX[,!colnames(BIN.MX)%in%as.character(Cp.v)],
                               maxCp=apply(X=BIN.MX[,as.character(Cp.v)], MARGIN=1, FUN=function(rw){
                                 max(rw*Cp.v)
                               }), stringsAsFactors=FALSE
    )
    return(BIN.MX)
  })
  BIN.MX <- do.call("rbind", BIN.MX)
  rownames(BIN.MX) <- NULL
  #-------------------Feature overlap
  df <- list()
  for(i in 1:foi.v.len){
    
    ct <- ct.v[i]
    if(ct%in%Ct.v){
      ct.TF <- BIN.MX[[ct]]>0
    } else if(ct=="hg19"){
      # FC as go to ct
      ct.TF <- BIN.MX[["FC"]]>0
      print("Choosing FC as default ct for ct=hg19.", quote=FALSE)
    } else {
      stop("Invalid ct.", quote=FALSE)
    }
    
    foi <- foi.v[i]
    foi.bed <- read.table(file=paste0(foi.dir, "/", foi), stringsAsFactors=FALSE, 
                          header=FALSE)[,1:3]
    foi <- gsub(x=foi, pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
    foi.bed <- foi.bed[foi.bed[,1]%in%chr.v,]
    foi.len <- sum(foi.bed[,3]-foi.bed[,2]+1L)
    
    common.chr <- intersect(foi.bed[,1], BIN.MX$chr)
    fr.v <- sapply(X=c(0, Cp.v), simplify=TRUE, FUN=function(Cp){
      # intersect() issues a warning when there are seqnames not common to both ranges
      incl.TF <- BIN.MX[,"maxCp"]==Cp & ct.TF & BIN.MX$chr%in%common.chr
      if( sum(incl.TF)==0 ){
        return(0)
      } else {
        B <- makeGRangesFromDataFrame(df=BIN.MX[incl.TF,c("chr", "start", "end")], 
                                      seqnames.field="chr", start.field="start", 
                                      end.field="end"); rm(incl.TF)
      }
      A=makeGRangesFromDataFrame(df=foi.bed[foi.bed[,1]%in%as.character(seqnames(B)),], 
                                 seqnames.field="V1", start.field="V2", 
                                 end.field="V3")
      return( sum(width(GenomicRanges::intersect(A,B)))/foi.len )
    }); rm(foi.bed, ct.TF, foi.len, common.chr, ct); gc()
    names(fr.v) <- c(0,Cp.v)
    df[[foi]] <- cbind.data.frame( foi=foi, fr=fr.v, Cp=as.character(c(0, Cp.v)) )
    rm(fr.v)
    print(paste0(foi, " done!"), quote=FALSE)
    
  } # foi.v.len for loop end
  df <- do.call("rbind.data.frame", df)
  df$Cp <- factor(df$Cp, levels=as.character(c(0, Cp.v)))
  save(df, file=paste0(out.dir, "/", out.name, "_foicentric.RData"))
  
} else {
  load(file=paste0(out.dir, "/", out.name, "_foicentric.RData"))
}


coul <- colorRampPalette(rev(brewer.pal(n=11,name="Spectral")))(length(Cp.v))
coul <- c("gray50", coul)
names(coul) <- c(0, Cp.v)
#-------------------Plot
ggplot(data=df, aes(x=foi, y=fr, fill=Cp)) +
  geom_col(position="fill") +
  scale_fill_manual(values=coul[levels(df$Cp)]) +
  labs( title=paste0(out.name), x="", y="Fraction of feature total bp",
        fill=expression(bold("c"["p"])) ) + 
  bgr2 +
  theme(axis.text.x=element_text(angle=90, size=5, hjust=1))
 
ggsave(filename=paste0(out.dir, "/", out.name, "_foicentric.pdf"),
       height=20, width=50, units="in", limitsize=FALSE)

# rm(list=ls()); gc()
