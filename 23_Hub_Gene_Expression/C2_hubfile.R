################################################################################
# Make hub files (similar to downloaded output from shiny app) for each hub 
# center provided
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(wk.dir, "/out_basePersist")
#persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_hubfile")
hubcenterPath = paste0(out.dir, "/All_hubcenter.csv")
annofilePath = paste0(wk.dir, "/txTable/hg19annoLTr_ALL")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v =paste("chr", c(1:22, "X"), sep="")
bin.len = 40000
topCP = 3; cp.v = 1:21
# If ct = "All", no filtering based on cell/tissue.
ct = "All"; ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
                    "AG", "Ov", "Bl", "MesC","MSC", "NPC", "TLC", "ESC", "FC", "LC")
# Gap in terms of bin
gap.bin = 50
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(GenomicRanges)
source(paste0(wk.dir, "/lib/GEN_WhichOverlap.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
ct.id <- NULL
if( !is.null(ct) ){ ct.id <- ct; ct.id <- paste0(ct.id[!is.na(ct.id)], "_") }
out.id <- paste0(ct.id, "topCP", topCP, "_gapBin", gap.bin)

anno.df <- read.delim(file=annofilePath, stringsAsFactors=FALSE, header=TRUE)
anno.df <- data.frame(chr=anno.df$chrom, txStart=anno.df$txStart, txEnd=anno.df$txEnd,
                      gene=anno.df$name2, accession=anno.df$name, uniqueID=anno.df$uniqueID,
                      stringsAsFactors=FALSE)

hcenter <- read.csv(file=hubcenterPath, header=TRUE, stringsAsFactors=FALSE,
                    check.names=FALSE, na.strings="")
hcenter <- hcenter[,c("chr", as.character(gap.bin))]

for(chr in chr.v){
  
  center.v <- hcenter[hcenter$chr==chr, as.character(gap.bin)]
  if( is.na(center.v) ){
    print(paste0("No center for ", chr, " gapBin=", gap.bin), quote=FALSE)
    next
  }
  center.v <- as.numeric( unique(strsplit(x=center.v, split=";")[[1]]) )
  center.v.len <- length(center.v)
  
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, "_topCP4.RData"))
  #load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
  incl.TF <- PERSIST.MX$ntis%in%tail(x=cp.v, n=topCP) & (PERSIST.MX$hits$j-PERSIST.MX$hits$i) > gap.bin
  if(ct%in%ct.v){
    incl.TF <- incl.TF & PERSIST.MX$hits[[ct]]>0
    print("Filtering based on cell/tissue.", quote=FALSE)
  } else if(ct=="All"){
    print("No filtering based on cell/tissue.", quote=FALSE)
  } else {
    stop("Invalid ct argument.")
  }
  ij.df <- PERSIST.MX$hits[incl.TF,c("i", "j")]
  rm(incl.TF, PERSIST.MX)
  
  anno.TF <- anno.df$chr==chr
  for(cl in 1:center.v.len){
    center <- center.v[cl]
    ij.TF <- ij.df$i==center | ij.df$j==center
    mem.sum <- table(c(ij.df$i[ij.TF], ij.df$j[ij.TF]))
    rm(center, ij.TF)
    mem.v <- as.character( sort(as.numeric(names(mem.sum)), decreasing=FALSE) )
    #-------------------Identify genes overlapping
    bin.end <- as.numeric(mem.v)*bin.len
    hits.mx <- WhichOverlap(start.query=bin.end-bin.len+1, 
                            end.query=bin.end, 
                            space.query=rep(x=chr, times=length(mem.v)),
                            start.subject=anno.df$txStart[anno.TF],
                            end.subject=anno.df$txEnd[anno.TF],
                            space.subject=anno.df$chr[anno.TF], 
                            maxgap=-1L, minoverlap=1L, type="any")
    rm(bin.end)
    if(nrow(hits.mx)==0){
      print(paste0("No overlapping genes for hub ", cl, "."), quote=FALSE)
      temp <- data.frame(gene=NA, accession=NA, uniqueID=NA)
    } else {
      temp <- by(INDICES=mem.v[hits.mx[,"query"]],
                 data=anno.df[anno.TF,c("gene", "accession", "uniqueID")][hits.mx[,"subject"],],
                 FUN=function(mx){
                   apply(X=mx, MARGIN=2, FUN=paste, collapse=";")
                 })
      temp <- cbind.data.frame(bin=names(temp), do.call(rbind, temp), stringsAsFactors=FALSE)
      temp <- temp[match(mem.v, temp$bin), -1]
    }
    
    df <- data.frame(bin=mem.v, temp, partner=NA, Nij=as.numeric(mem.sum[mem.v]), 
                     stringsAsFactors=FALSE)
    rm(temp, mem.sum, mem.v)
    write.csv(df, file=paste0(out.dir, "/", gcb, "_", chr, "_", out.id, "_hub", cl, ".csv"), 
              quote=FALSE, row.names=FALSE)
    
    # Save genelist for functional annotation
    genelist <- unlist(strsplit(x=df$gene, split=";"))
    genelist <- unique(genelist[!is.na(genelist)])
    writeLines(text=genelist, 
               con=paste0(out.dir, "/genes_txt/", gcb, "_", chr, "_", out.id, "_hub", cl, "_genes.txt"))
    
    rm(df); gc()
    print(paste0("hub ", cl), quote=FALSE)
  }
  rm(ij.df, anno.TF); gc()
  
  print(paste0(chr, " gapBin=", gap.bin, " done!"), quote=FALSE)
  
}

# rm(list=ls()); gc()



