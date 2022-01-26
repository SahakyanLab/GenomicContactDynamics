################################################################################
# Make plot showing the length of contact bins after hg19 -> hg38 liftover
# and the percentage of bins per chromosome that were split into multiple ranges
# after liftover.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/6_Location_Chrom3D"
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
model.id    = "H1-hESC_LMNB1_hg38"
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
LOdata.dir  = paste0(wk.dir, "/out_mapToHiCcontactBins_Chrom3D/H1-hESC_LMNB1_hg38") 
out.dir     = paste0(wk.dir, "/out_binsAfterLO_plot") 
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt") 
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
bin.len = 40000
chr.v = paste0("chr", c(1:22, "X"))
LOchain = "hg19ToHg38"
plotOnly = T
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste(gcb, LOchain, model.id, bin.len, sep="_")

if(plotOnly==F){
  
  chrlen.df <- read.delim(file=chrlen.file, header=T, sep="\t", stringsAsFactors=F)
  totbin.chr.v <- setNames(object=chrlen.df$bins.40kb, nm=chrlen.df$chromosome)
  
  if( !identical(as.numeric(totbin.chr.v), as.numeric(chrlen.df$bins.40kb)) ){
    stop("Checkpoint 1.")
  }
  rm(chrlen.df)
  
  totbin.ij.v <- setNames(object=rep(NA, times=length(chr.v)), nm=chr.v)
  totbin.1range.kept.len.LO.v <- totbin.ij.LO.v <- totbin.ij.v
  LO <- list()
  for(chr in chr.v){
    
    # Get total contact bins per chr
    load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    totbin.ij.v[[chr]] <- length(unique( c(PERSIST.MX$hits$i, PERSIST.MX$hits$j) ))
    
    rm(PERSIST.MX)
    gc()
    
    # Contains contact bins lifted over
    load(paste0(LOdata.dir, "/", chr, "_Persist_min2Mb_", LOchain, "_RAW.RData"))
    totbin.ij.LO.v[[chr]] <- length(unique(LO.mx$group))
    
    # This can happen if a contact bin was not converted at all
    if( totbin.ij.v[[chr]]!=totbin.ij.LO.v[[chr]] ){
      warning(paste0(chr, ": Total contact bins converted not equal to total contact bins."))
    }
    
    # Remove bins split into multiple ranges after lift over
    LO.mx <- LO.mx[!(duplicated(x=LO.mx$group, fromLast=F) | duplicated(x=LO.mx$group, fromLast=T)),
                   c("group", "seqnames", "width")]
    
    if( any(duplicated(LO.mx$group)) | any(LO.mx$group>as.numeric(totbin.chr.v[chr])) ){
      stop(paste0(chr, ": Checkpoint 2."))
    }
    
    if( any(LO.mx$width>bin.len) ){
      stop(paste0(chr, ": Length greater than bin.len."))
    } else {
      totbin.1range.kept.len.LO.v[[chr]] <- sum(LO.mx$width==bin.len)
    }
    
    LO[[chr]] <- LO.mx
    
    rm(LO.mx)
    
    print(paste0(chr, " done!"), quote=F)
    
  }
  
  LO <- do.call("rbind.data.frame", LO)
  rownames(LO) <- NULL
  
  LO$seqnames <- as.character(LO$seqnames)
  LO$group <- as.numeric(as.character(LO$group))
  LO$width <- as.numeric(as.character(LO$width))
  
  print(range(LO$width, na.rm=F))
  #[1]   324 40000
  
  if( any(!chr.v%in%LO$seqnames) | any(!LO$seqnames%in%chr.v) ){
    stop("Missing chromosomes.")
  }
  
  save(LO, file=paste0(out.dir, "/", out.name, "_widthDensAfterLO.RData"))
  
} else {
  load(file=paste0(out.dir, "/", out.name, "_widthDensAfterLO.RData"))
}

# Out of valid bins, plot percentage of those equal to and less than the bin.len

num.v = c(sum(LO$width==40000), sum(LO$width<40000))
perc.v = setNames(object=num.v/sum(num.v)*100, 
                  nm=c("= 40 kb", "< 40 kb"))

pdf(file=paste0(out.dir, "/", out.name, "_widthDensAfterLO_pie.pdf"),
    height=10, width=20)
par(mfrow=c(1,2))

pie(x=num.v, labels=names(perc.v), col=c("gray90", "#55bde6"), cex.main=0.5,
    main=paste0(out.name, "afterLO \n", 
                paste(paste(names(perc.v), perc.v, sep="="), collapse="; "))
    )

plot(hist(x=LO$width, breaks=seq(0,40000, 5000), col="#55bde6"))

dev.off()

# Store numbers

totbin.1range.LO.v <- table(LO$seqnames)
out.df <- data.frame(chr=as.character(chr.v), 
                     totbin.chr=as.numeric(totbin.chr.v[chr.v]), 
                     totbin.ij=as.numeric(totbin.ij.v[chr.v]), 
                     totbin.ij.LO=as.numeric(totbin.ij.LO.v[chr.v]),
                     totbin.1range.LO=as.numeric(totbin.1range.LO.v[chr.v]), 
                     totbin.1range.kept.len.LO=as.numeric(totbin.1range.kept.len.LO.v[chr.v]),
                     stringsAsFactors=F)
out.df <- rbind.data.frame(out.df, c(chr=1, colSums(out.df[,-1])),
                           stringsAsFactors=F)
out.df$chr[nrow(out.df)] <- "chrALL"

tmp <- as.numeric(apply(X=data.matrix(out.df[,-1]), MARGIN=1, FUN=diff))
if( any(tmp>0) ){
  stop("Numbers not decreasing as you go from right to left columns.")
} else {
  
  write.csv(x=out.df, row.names=F,
            file=paste0(out.dir, "/", out.name, "_countofvalidbinsafterLO.csv"))
  
}

# rm(list=ls()); gc()

# Out of valid bins, plot percentage of those equal to and less than the bin.len
#pdf(file=paste0(out.dir, "/", out.name, "_widthDensAfterLO_pie.pdf"),
#    height=10, width=10)

#pie.val <- table(LO$width==bin.len)
#pie.val <- pie.val/sum(pie.val)
#pie.val <- pie.val[c("TRUE", "FALSE")]
#me <- paste(c("bin.len\n", "less than bin.len\n"), 
#             paste0(format(x=pie.val*100, digits=2), " %"))
#names(pie.val) <- nme
#pie(pie.val, col=rev(c("deepskyblue3", "gray80")))

#dev.off()
