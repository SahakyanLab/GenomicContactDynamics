################################################################################
# Compare discordance of original and shuffled HiC contacts; contacts can be
# selected based on overlap with a feature
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",

# Avoid left to right partial matching by $
options(warnPartialMatchDollar = T)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    lib = paste0(home.dir, "/DPhil/lib")
    data.dir = paste0(home.dir, "/Database")
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/12_Shuffling")
    orig.dir = paste0(home.dir, "/DPhil/GCD_polished/11_Complementarity/z_ignore_git/out_constraints/merged_final")
    shuff.dir = paste0(wk.dir, "/z_ignore_git/out_constraints/merged_final")
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/t1-data/user/ltamon"
    lib = paste0(home.dir, "/DPhil/lib")
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/8_ShuffleContactBins")
    data.dir = paste0(home.dir, "/Database")
    orig.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/pending/11_Constraints/out_constraints/merged_final")
    shuff.dir = paste0(wk.dir, "/out_constraints/merged_final")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
#persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/persist_HiCNorm")
out.dir = paste0(wk.dir, "/out_compare_new")
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/raw")
foifile = paste0(wk.dir, "/foifile/foifile6")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(1:22, "X"), sep="")
gcb = "min2Mb"
bin.len = 40000L
kmer.len = 7L
type = "align" # kmer | align
affix = "_ijShuffled"
plotOnly = F
filterByFoi = T
filterByCelltype = T
mannwhit = T
out.id = "chrALL_6"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(compiler)
library(data.table)
library(GenomicRanges)
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
### FUNCTION ###################################################################
HicHybridCompare <- function(orig.dir = '/dir',
                             shuff.dir = '/dir',
                             out.dir = '/dir',
                             foi.dir = 'feature directory',
                             foifile = 'feature name',
                             gcb = 'contact gap id',
                             chr.v = 'chromosome/s whose data will be combined',
                             bin.len = "contact bin resolution",
                             tot.len.v = "total chr length; same order as chr.v",
                             kmer.len = "for c||kmer",
                             # Identifier of the shuffled set
                             affix = "_ijShuffled",
                             plotOnly = 'T/F',
                             filterByFoi = 'T/F',
                             mannwhit = 'T/F',
                             out.id = '<chr>_<affix>'
                             ){
  
  out.id <- paste0(out.id, "_filterByFoi", filterByFoi, "_filterByCelltype", 
                   filterByCelltype)
  
  if(filterByFoi){
    
    foi <- readLines(con=foifile)
    foi.df <- read.delim(file=paste0(foi.dir, "/", foi), sep="\t", stringsAsFactors=F, 
                         header=F)
    
    out.id <- paste0(out.id, "_", foi)
    ct <- strsplit(x=foi, split="_", fixed=T)[[1]][2]
    
  }
  
  if(plotOnly==F){
    
    HYBCOMB.DF <- list()
    
    for(chr in chr.v){
      
      load(file=paste0(orig.dir, "/", chr, "_", type, "_", gcb, ".RData"))
      CII.MX <- CII.MX[!is.na(CII.MX[,"Cp"]),]
      
      # Filter by cell type if required
      if(filterByCelltype & filterByFoi){
        
        print(paste0(out.id, ": Account cell type of ", foi, " feature."), quote=F)
        
        load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
        ct.v <- colnames(PERSIST.MX$hits)[!colnames(PERSIST.MX$hits)%in%c("i", "j")]
        
        if(ct%in%ct.v){
          
          if( !identical(as.numeric(CII.MX[,"i"]), as.numeric(PERSIST.MX$hits$i)) |
              !identical(as.numeric(CII.MX[,"j"]), as.numeric(PERSIST.MX$hits$j)) ){
            stop( paste0(out.id, ": CII.MX and PERSIST.MX$hits rownames don't match.") )
          }
          
          CII.MX <- CII.MX[PERSIST.MX$hits[[ct]]>0,]
          
        } # ct.v
        
        rm(PERSIST.MX, ct.v)
        
      } # filterByCelltype
      
      rownames(CII.MX) <- NULL
      HYBCOMB.DF[[chr]] <- data.frame(chr=chr, CII.MX, set="1", stringsAsFactors=F)
      rm(CII.MX)
      
      # Shuffled
      load(file=paste0(shuff.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"))
      rownames(CII.MX) <- NULL
      CII.MX <- CII.MX[!is.na(CII.MX[,"Cp"]),]
      HYBCOMB.DF[[chr]] <- rbind(HYBCOMB.DF[[chr]], data.frame(chr=chr, CII.MX, set="2",
                                                               stringsAsFactors=F))
      
    }
    
    HYBCOMB.DF <- do.call("rbind", HYBCOMB.DF)
    rownames(HYBCOMB.DF) <- NULL
    data.table::setnames(x=HYBCOMB.DF, old="C..", new="C||", skip_absent=F)
    save(HYBCOMB.DF, file=paste0(out.dir, "/", out.id, "_HybComb_", type, "_", gcb, ".RData"))
    
  } else {
    load(paste0(out.dir, "/", out.id, "_HybComb_", type, "_", gcb, ".RData"))
  }
  
  #-------------------If required, select only contacts in HYBCOMB.DF overlapping with feature
  
  if(filterByFoi){
    
    print(paste0(out.id, ": Filtering contacts by ", foi, " feature."), quote=F)
    
    names(tot.len.v) <- chr.v
    
    for(chr in chr.v){
      
      fchr.TF <- foi.df$V1==chr
      
      if( sum(fchr.TF)==0 ){
        next
      } else {
        
        tot.bin <- ceiling(tot.len.v[[chr]]/bin.len)
        bin.end <- (1:tot.bin)*bin.len
        bin.start <- bin.end-bin.len+1L
        bin.end[tot.bin] <- tot.len.v[[chr]]
        
        olap <- WhichOverlap(start.query=bin.start, 
                             end.query=bin.end, 
                             space.query=rep(x=chr, times=tot.bin),
                             start.subject=foi.df[fchr.TF,2], 
                             end.subject=foi.df[fchr.TF,3], 
                             space.subject=foi.df[fchr.TF,1],
                             maxgap=-1L, minoverlap=1L,
                             type="any")
        ubins.ind <- as.numeric(unique(olap[,1]))
        rm(olap, tot.bin, bin.end, bin.start)
        
        hchr.TF <- HYBCOMB.DF$chr==chr
        drop.TF <- hchr.TF & !(HYBCOMB.DF$i%in%ubins.ind & HYBCOMB.DF$j%in%ubins.ind)
        HYBCOMB.DF <- HYBCOMB.DF[!drop.TF,]
        
        rm(hchr.TF, drop.TF, ubins.ind)
        
      }
      
      rm(fchr.TF)
      
    } # chr.v for loop end
    
  }
  
  # Change C|| to CII to make it valid for the Mann-Whitney test
  data.table::setnames(x=HYBCOMB.DF, old="C||", new="CII", skip_absent=F)
  
  #-------------------Plotting requirements
  
  clnme <- colnames(HYBCOMB.DF)
  feat.v <- clnme[ !clnme%in%c("chr", "i", "j", "Cp", "set") ]
  
  # Plot labels
  if(type=="align"){
    ylab <- list( CII=bquote(bold( "c"["||"]~"align" )) )
  } else if(type=="kmer"){
    
    ylab <- list(CII=bquote(bold( "c"["||"]~"kmer" )),
                 Gfree=bquote(bold( "Gfree, kcal/mol" )),
                 sdDifference=bquote(bold( "s ("~"c"["||"]~")" ))
    )
    
  } else {
    stop("Invalid type.")
  }
  
  id <- paste0(out.id, "_" , gcb, "_", type)
  
  # Specify order to make sure that orig will always be on the left and shuff on the right
  HYBCOMB.DF$set <- factor(x=as.character(HYBCOMB.DF$set), levels=c("1", "2"))
  
  # Turn into a factor, whose levels are ordered from 1 to 21. This order determines the
  # sequence along the x-axis.
  cp.v <- sort( as.numeric(as.character(unique(HYBCOMB.DF$Cp))), decreasing=F )
  cp.v.len <- length(cp.v)
  HYBCOMB.DF$Cp <- factor(as.character(HYBCOMB.DF$Cp), levels=as.character(cp.v))
  
  # Initialize output matrix with the p-values
  if(mannwhit){
    pval.mx <- matrix(data=NA, nrow=cp.v.len, ncol=length(feat.v), 
                      dimnames=list(cp.v, feat.v))
  }
  
  #-------------------Record important numbers
  
  # Count non-NA contacts orig and shuff to check the difference
  if( any(is.na(HYBCOMB.DF$Cp)) ){
    stop(paste0(out.id, ": Checkpoint."))
  }
  ij.orig.TF <- as.character(HYBCOMB.DF$set)=="1"
  ij.shuff.TF <- !ij.orig.TF
  feat.nonNA.TF <- !is.na(HYBCOMB.DF$CII)
  
  # Number of contacts
  tot.ij.orig <- sum(ij.orig.TF)
  perc.ij.shuff <- sum(ij.shuff.TF)/tot.ij.orig*100
  perc.nonNA.ij.orig <- sum(ij.orig.TF & feat.nonNA.TF)/tot.ij.orig*100
  perc.nonNA.ij.shuff <- sum(ij.shuff.TF & feat.nonNA.TF)/tot.ij.orig*100
  
  perCp.orig <- table(HYBCOMB.DF$Cp[feat.nonNA.TF & ij.orig.TF])
  rm(ij.orig.TF)
  #percPerCp.orig <- perCp.orig/tot.ij.orig*100
  percPerCp.shuff <- (table(HYBCOMB.DF$Cp[feat.nonNA.TF & ij.shuff.TF])/perCp.orig)*100
  rm(ij.shuff.TF, feat.nonNA.TF)
  
  mx <- cbind(orig=c(tot=tot.ij.orig, nonNA=perc.nonNA.ij.orig, perCp.orig),
              shuff=c(tot=perc.ij.shuff, nonNA=perc.nonNA.ij.shuff, percPerCp.shuff))
  write.table(mx, file=paste0(out.dir, "/", id, "_counts"), col.names=T, row.names=T,
              quote=F, sep="\t")
  
  rm(perCp.orig, percPerCp.shuff)
  gc()
  
  #-------------------Plot
  
  for(feat in feat.v){
    
    jpeg(filename=paste0(out.dir, "/", id, "_", feat, "_combBP.jpeg"),
         units="in", width=12, height=10, res=500)
    # default = par(mar=c(5,4,4,2)+0.1);  c(bottom, left, top, right)
    # default = par(mgp=c(3, 1, 0)); c(axis title, axis labels, axis line)
    par(mar=c(5.5, 6.5, 6.5, 2)+0.1, mgp=c(3, 1.5, 0)) 
    myplot <- boxplot(as.formula(paste0(feat,"~set*HYBCOMB.DF$Cp")), outline=F,
                      data=HYBCOMB.DF, boxwex=0.6, xlab="", ylab="", main="", 
                      cex.axis=2.5, col=c("#FDC776" , "gray91"), xaxt="n")
    # Labelling x-axis
    # "1" "1" "2" "2" "3" "3"
    xlabs <- sapply( strsplit(x=myplot$names , split='\\.') , 
                     FUN=function(x) x[[2]] )
    # "1" "2" "3" "4"
    xlabs <- xlabs[seq(1, length(xlabs), 2)]
    axis(side=1, at=seq(1.5, cp.v.len*2 , by=2), 
         labels=1:21, tick=T, cex.axis=2.5)
    
    # X-axis title
    mtext(side=1, text=expression( bold("c"["p"]) ), line=4.5, cex=3)
    # Y-axis title
    mtext(side=2, text=bquote(bold( .(ylab[[feat]]) )), line=3.5, cex=3)
    # Plot title
    mtext(side=3, text=paste0(id, "_outline=NA"), line=1.5, cex=0.5)
    
    # Line separating two boxplots per Cp
    for(i in seq(0.5 , cp.v.len*2+1, 2)){ abline(v=i, lty=1, col="gray50")}
    legend("topright", legend=c( expression(bold("orig")), expression(bold("shuff")) ), 
           col=c("#FDC776" , "gray91"), pch=15, bty="o", bg="white", pt.cex=3, 
           cex=2, horiz=F, inset=c(0,-0.17), xpd=T)
    dev.off()
    
    #---------------------------------------
    
    if(mannwhit){
      
      # Wilcoxon-Mann-Whitney test (non-parametric, unpaired/independent groups) 
      # to compare orig and shuff values per Cp
      # The test assumes that the shape of the two distributions are similar so
      # check the boxplots
      # If both x and y are given and paired is FALSE, a Wilcoxon rank sum test
      # (equivalent to the Mann-Whitney test) is carried out.
    
      for(cp in cp.v){
        
        cp <- as.character(cp)
        mw <- wilcox.test(as.formula(paste0(feat, "~set")), alternative="two.sided",
                          data=HYBCOMB.DF[as.character(HYBCOMB.DF$Cp)==cp & !is.na(HYBCOMB.DF[[feat]]), 
                                          c(feat, "set")],
                          paired=F)
        # tt <- t.test("set", feat)
        pval.mx[cp,colnames(pval.mx)==feat] <- mw$p.value
        rm(mw, cp)
        
      }
      
    }
   
    print(paste0(out.id, ": ", feat, " done!"), quote=F)
    
  } # feat.v for loop
  
  if(mannwhit){
    
    write.table(x=pval.mx, file=paste0(out.dir, "/", id, "_OlessS_mwtest"),
                col.names=T, row.names=T, quote=F, sep="\t")
    
  }
  
}
################################################################################
HicHybridCompare <- cmpfun(HicHybridCompare, options=list(suppressUndefined=T))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrlen.df <- read.delim(file=chrlen.file, header=T, stringsAsFactors=F,
                        row.names="chromosome")
tot.len.v <- setNames(object=chrlen.df[chr.v, "length.bp"], nm=chr.v)

HicHybridCompare(
  
  orig.dir=orig.dir,
  shuff.dir=shuff.dir,
  out.dir=out.dir,
  foi.dir=foi.dir,
  foifile=foifile,
  gcb=gcb,
  chr=chr.v,
  bin.len=bin.len,
  tot.len.v=tot.len.v,
  kmer.len=kmer.len, 
  affix=affix, 
  plotOnly=plotOnly,
  filterByFoi=filterByFoi,
  mannwhit=mannwhit,
  out.id=out.id
  
)

# rm(list=ls()); gc()


