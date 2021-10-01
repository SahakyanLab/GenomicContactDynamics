################################################################################
# Map feature to contacts, contact- and feature-centric wise.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    #wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/persist_HiCNorm")
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/raw")
foifile = paste0(wk.dir, "/foifile/foivsij/foifile26")
out.dir = paste0(wk.dir, "/out_foiVsij")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
bin.len = 40000
Cp.v = 1:21
Ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", 
         "SB", "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC",
         "FC", "LC")
out.id = "26" 
type.olap = "any" # "any" | "within"
query = "bin" #"bin" | "foi"
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(GenomicRanges)
library(ggplot2)
library(reshape)
#library(ggsci)
library(RColorBrewer)
#library(yarrr); base.pal <- yarrr::piratepal("basel")
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
source(paste0(lib, "/finaliseFOI.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, " done..."), quote=FALSE)
Cp.v <- sort(Cp.v, decreasing=FALSE)
#-------------------List of features and corresponding cell/tissue
foi.v <- unique(finaliseFOI(foi.dir=foi.dir, foifile=foifile))
foi.v.len <- length(foi.v)
foilab.v <- gsub(x=foi.v, pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
ct <- unlist(
  lapply(X=strsplit(x=foi.v, split="ct_|\\_"), FUN=function(x)x[2])
)
ct <- unique(ct)
# Features should be from the same ct
if( length(ct)!=1 ){ stop("Features from different cell/tissue.") }

out.name <- paste0(gcb, "_", ct, "_", out.id, "_", type.olap, "_query", query)

if(plotOnly==FALSE){
  
  #-------------------Prepare foi.bed
  foi.bed <- sapply(X=1:foi.v.len, simplify=FALSE, FUN=function(i){
    bed <- read.table(file=paste0(foi.dir, "/", foi.v[i]), stringsAsFactors=FALSE, 
                      header=FALSE)[,1:3]
    cbind.data.frame(bed, foi=foilab.v[i], stringsAsFactors=FALSE)
  })
  foi.bed <- do.call("rbind.data.frame", foi.bed)
  chr.v <- intersect(unique(foi.bed[,1]), chr.v)
  rownames(foi.bed) <- NULL
  #-------------------Prepare foi pairings as contact groups
  group.v <- 10^(1:(length(foilab.v)+2L))
  names(group.v) <- c("none", "promisc", foilab.v)
  
  pair.v <- group.v + group.v
  names(pair.v) <- paste(names(group.v), names(group.v), sep="-")
  comb <- combn(x=names(group.v), m=2, FUN=paste, collapse="-", simplify=FALSE)
  temp <- lapply(X=comb, FUN=function(pair){
    sum(group.v[ strsplit(x=pair, split="-")[[1]] ])
  })
  names(temp) <- comb; pair.v <- c(pair.v, temp); rm(temp)
  if( any(duplicated(pair.v)) ){ stop("Duplicate in pair.v.") }
  #-------------------Output matrix
  FEATVSCP.MX <- matrix(data=0, nrow=length(Cp.v), ncol=length(pair.v)+1L, 
                        dimnames=list(Cp.v, c("All", pair.v))
  )
  #-------------------Gather data per
  for(chr in chr.v){
    
    chr.TF <- foi.bed[,1]==chr
    
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    if(ct%in%Ct.v){
      ct.TF <- PERSIST.MX$hits[[ct]]>0
    } else if(ct=="hg19"){
      ct.TF <- rep(TRUE, times=nrow(PERSIST.MX$hits))
    } else {
      stop("Invalid ct.")
    }
    ubins <- unique(c( unique(PERSIST.MX$hits$i[ct.TF]), unique(PERSIST.MX$hits$j[ct.TF]) ))
    ubins <- sort(ubins, decreasing=FALSE)
    ij.df <- cbind(PERSIST.MX$hits[ct.TF, c("i", "j")], Cp=PERSIST.MX$ntis[ct.TF])
    rm(ct.TF, PERSIST.MX); gc()
    
    # Overlap
    ubin.end <- ubins*bin.len
    ubin.start <- ubin.end-bin.len+1L
    
    if( type.olap=="any" | (type.olap=="within" & query=="bin") ){
      
      olap.df <- WhichOverlap(start.query=ubin.start, 
                              end.query=ubin.end, 
                              space.query=chr,
                              start.subject=foi.bed[chr.TF,2], 
                              end.subject=foi.bed[chr.TF,3], 
                              space.subject=foi.bed[chr.TF,1],
                              maxgap=-1L, minoverlap=1L,
                              type=type.olap)
      print("query-bin; subject-foi", quote=FALSE)
      
    } else if( type.olap=="within" & query=="foi"){
      
      olap.df <- WhichOverlap(start.subject=ubin.start, 
                              end.subject=ubin.end, 
                              space.subject=chr,
                              start.query=foi.bed[chr.TF,2], 
                              end.query=foi.bed[chr.TF,3], 
                              space.query=foi.bed[chr.TF,1],
                              maxgap=-1L, minoverlap=1L,
                              type=type.olap)
      print("subject-bin; query-foi", quote=FALSE)
      # Switch to query-bin, subject-foi format
      colnames(olap.df) <- c("subject", "query")
      
    } else {
      stop("Invalid type.olap and query arguments.")
    }
    
    olap.df <- as.data.frame(olap.df)
    rm(ubin.end, ubin.start)
    
    # Classify unique bins (FIX DUPLICATED BIN)
    olap.df$subject <- foi.bed[chr.TF, "foi"][olap.df$subject]
    olap.df$query <- ubins[olap.df$query]
    temp <- by(data=olap.df$subject, INDICES=olap.df$query, simplify=FALSE,
               FUN=function(x)unique(as.character(x)))
    promisc.bin <- as.numeric( names(temp)[lapply(X=temp, FUN=length)>1] ); rm(temp)
    olap.df <- olap.df[!olap.df$query%in%promisc.bin,] 
    olap.df <- olap.df[!duplicated(olap.df$query),]
  
    ubin.grp <- rep(NA, times=length(ubins))
    names(ubin.grp) <- ubins
    ubin.grp[ as.character(olap.df$query) ] <- olap.df$subject
    ubin.grp[ as.character(promisc.bin) ] <- "promisc"
    ubin.grp[ is.na(ubin.grp) ] <- "none"
    if( any(is.na(ubin.grp)) ){ stop("Unclassified bin.") }
    
    # Classify contacts
    ij.df$group <- group.v[ ubin.grp[as.character(ij.df$i)] ] + group.v[ ubin.grp[as.character(ij.df$j)] ]
    if( any(is.na(ij.df$group)) ){ stop("Unclassified contact.") }
    
    for(Cp in Cp.v){
      cp.TF <- ij.df$Cp==Cp; Cp <- as.character(Cp)
      v <- table(ij.df$group[cp.TF])
      FEATVSCP.MX[Cp, c("All", names(v))] <- FEATVSCP.MX[Cp, c("All", names(v))] + c(sum(cp.TF), unname(v))
      rm(v, Cp)
    }
    print(paste0(chr, " done!"), quote=FALSE)
    
  } # chr.v for loop end
  colnames(FEATVSCP.MX) <- c("All", names(pair.v))
  save(FEATVSCP.MX, file=paste0(out.dir, "/", out.name, "_foiVsij.RData"))
  
} else {
 load(file=paste0(out.dir, "/", out.name, "_foiVsij.RData")) 
}
v <- colSums(FEATVSCP.MX)
if( v["All"]!=sum(v[-1]) ){ stop("Checkpoint.") } 
ij.num <- paste(paste(names(v), v, sep=":"), collapse="_")
#-------------------Contact-centric plot
df <- melt.array(FEATVSCP.MX[,-1]/FEATVSCP.MX[,"All"])
colnames(df) <- c("Cp", "Pair", "Value")
df$Pair <- factor(df$Pair, levels=unique(df$Pair))
x.lab <- Cp.v; x.lab[(x.lab%%2)==0] <- ""

ggplot(data=df, aes(x=Cp, y=Value, fill=Pair)) +
  geom_col(position="fill") +
  scale_x_continuous(breaks=Cp.v, labels=x.lab) + 
  #scale_fill_npg() + 
  labs(title=paste0(out.name, "_olaptype=", type.olap), 
       x=expression(bold("c"["p"])),
       y="Fraction of contacts", fill=NULL) + 
  bgr2 +
  theme(legend.text=element_text(size=5, face="bold")) 
ggsave(filename=paste0(out.dir, "/", out.name, "_foiVsij_ij.pdf"),
       width=10, height=10)
#-------------------Feature-centric plot
totij.pair <- colSums(FEATVSCP.MX)
df <- apply(X=FEATVSCP.MX, MARGIN=1, FUN=function(rw){
  v <- rw/totij.pair; v[is.nan(v)] <- 0
  return(v)
}); rm(FEATVSCP.MX)
totij.pair[totij.pair>0] <- 1
if( !identical(rowSums(df), totij.pair) ){ stop("Fractions don't add up to 1.") }

df <- melt.array(df)
colnames(df) <- c("Pair", "Cp", "Value")
df$Cp <- factor(df$Cp, levels=as.character(Cp.v))
df <- df[df$Value!=0,]
df$Pair <- gsub(x=df$Pair, pattern="-", replacement="\n", fixed=TRUE)
df$Pair <- factor( df$Pair, levels=unique(c("All", as.character(unique(df$Pair)))) )
coul <- colorRampPalette(rev(brewer.pal(n=11,name="Spectral")))(length(Cp.v))

ggplot(data=df, aes(x=Pair, y=Value, fill=Cp)) +
  geom_col(position="fill") +
  scale_fill_manual(values=coul) +
  labs(title=paste0(out.name, "_olaptype=", type.olap, "\n",
                    ij.num), 
       x="", y="Fraction of contacts", 
       fill=expression(bold("c"["p"]))
       ) + 
  bgr2 +
  theme(axis.text.x=element_text(angle=90, size=10, hjust=1),
        plot.title=element_text(size=2))
ggsave(filename=paste0(out.dir, "/", out.name, "_foiVsij_foi.pdf"),
       width=10, height=10)

# rm(list=ls()); gc()
