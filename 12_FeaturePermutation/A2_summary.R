################################################################################
# Summarise results of association in a plot of p-values per feature
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
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
    binmx.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/19_Circos/out_bindata"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
permt.dir = paste0(wk.dir,"/out_associate_runA_cp21_b1b2b3")
out.dir = paste0(wk.dir,"/out_summary_runA_cp21_b1b2b3")
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/reduced_b1b2b3")
foifile = paste0(wk.dir, "/foifile/foifile_b1b2b3")
### OTHER SETTINGS #############################################################
out.name = "runA_cp21_b1b2b3" #"ABscomp"
gcb = "min2Mb"
# Restrict chromosomes to:
chr.v = paste("chr", c(1:22, "X"), sep="")
Cp.v = 21
# Cell/tissue of features
ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
         "AG", "Ov", "Bl", "MesC","MSC", "NPC", "TLC", "ESC", "FC", "LC") #c("hg19", "FC", "ESC", "LC")
# Permutation
NTIMES = 10000
SEED = 342
# Selected evaluation functions
eval.f.v = c(`% number of A regions overlapping`="numOlapA",
             #`% number of A regions overlapping within`="numOlapAwithin",
             `% number of B regions overlapping`="numOlapB",
             #`% number of B regions overlapping within`="numOlapBwithin",
             `total length of intersection`="comOlap",
             `% intersection relative to A`="comOlapA",
             `% intersection relative to B`="comOlapB",
             `mean distance A from B`="meandist"
             #,
             #`mean length of A regions overlapping`="meanLenOlapA"
)
#col.v <- adjustcolor(c("#3F007D", "#6A51A3", "#9E9AC8", "gray90"), alpha.f=0.7)
col.v = c(greater="#6A51A3B3", less="#CCCCCCB3")
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/finaliseFOI.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# For calculating length of contacting regions (cell-type specific); 
# to be divided to comOlap
BIN.MX <- sapply(X=chr.v, simplify=FALSE, 
                 FUN=function(chr){
  load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
  BIN.MX <- cbind.data.frame(chr=chr, 
                             BIN.MX[,colnames(BIN.MX)%in%c("start", "end", ct.v, Cp.v)], 
                             stringsAsFactors=FALSE)
  # Select regions with Cp in Cp.v
  BIN.MX <- BIN.MX[rowSums(as.matrix(BIN.MX[,as.character(Cp.v)]))>0,]
})
BIN.MX <- do.call("rbind", BIN.MX)
rownames(BIN.MX) <- NULL
#totlen.bin0 <- sum(BIN.MX$end-BIN.MX$start+1)
#---------------------------------------
evalf.lab <- names(eval.f.v)
evalf.lab.len <- length(evalf.lab)
evalf.nme <- eval.f.v
out.id <- paste0("nperm", NTIMES, "_cp", Cp.v, "_seed", SEED)
#---------------------------------------
if(plotOnly==FALSE){
  
  foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)
  foi.v.len <- length(foi.v)
  
  # Initialize output tables
  PERMTSUM <- list()
  PERMTSUM$pval <- matrix(data=NA, ncol=length(eval.f.v), nrow=foi.v.len,  
                          dimnames=list(gsub(x=foi.v, pattern="ct\\_|foi\\_|desc\\_|\\.bed", 
                                             replacement=""), 
                                        evalf.lab))
  PERMTSUM$obs <- PERMTSUM$pval
  PERMTSUM$alt <- PERMTSUM$pval
  PERMTSUM$totlen <- matrix(data=NA, ncol=2, nrow=foi.v.len, 
                            dimnames=list(gsub(x=foi.v, pattern="ct\\_|foi\\_|desc\\_|\\.bed", 
                                               replacement=""), 
                                          c("foi", "bin")))
  PERMTSUM$foifile <- foi.v
                         
  for(i in 1:foi.v.len){
    foi <- foi.v[i]
    foi.bed <- read.table(file=paste0(foi.dir, "/", foi), stringsAsFactors=FALSE, 
                          header=FALSE)[,1:3]
    foi.bed <- foi.bed[foi.bed[,1]%in%chr.v,]
    foichr.v <- unique(foi.bed[,1])
    # foi.bed is 1-based and reduced
    totlen.foi <- sum(foi.bed[,3]-foi.bed[,2]+1)
    rm(foi.bed); gc()
    
    foi <- gsub(x=foi, pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
    ct <- strsplit(x=foi, split="_", fixed=TRUE)[[1]][1]
    incl.TF <- BIN.MX$chr%in%foichr.v
    if(ct%in%ct.v){ incl.TF <- incl.TF & BIN.MX[[ct]]==1 } 
    totlen.bin <- sum(BIN.MX[incl.TF,"end"]-BIN.MX[incl.TF,"start"]+1) 
    rm(incl.TF)
    
    # Load permutation test object, PERMT 
    load(file=paste0(permt.dir, "/", gcb, "_", foi, "_", out.id, "_permtest.RData"))
    for(x in 1:evalf.lab.len){
      eval <- evalf.lab[x]
      PERMTSUM$pval[foi,eval] <- as.numeric( PERMT[[eval]][["pval"]] )
      PERMTSUM$obs[foi,eval]  <- as.numeric( PERMT[[eval]][["observed"]] )
      PERMTSUM$alt[foi,eval]  <- PERMT[[eval]][["alternative"]] 
      if(x==1){
        PERMTSUM$totlen[foi,"foi"] <- totlen.foi
        PERMTSUM$totlen[foi,"bin"] <- totlen.bin
        rm(totlen.bin, totlen.foi)
      }
    }
    rm(PERMT); gc()
    print(paste0(foi, " done!"), quote=FALSE)
  }
  save(PERMTSUM, file=paste0(out.dir, "/", out.id, "_", out.name, "_permtsum.RData"))
  
} else {
  # Load PERMTSUM
  load(file=paste0(out.dir, "/", out.id, "_", out.name, "_permtsum.RData"))
}

foi.v <- rownames(PERMTSUM$pval)
for(i in 1:evalf.lab.len){
  
  eval <- evalf.lab[i]
  df <- data.frame(foi=foi.v, obs=PERMTSUM$obs[,eval], alt=PERMTSUM$alt[,eval],
                   pval=PERMTSUM$pval[,eval], row.names=NULL, stringsAsFactors=FALSE)
  df$alt <- factor(df$alt, levels=c("greater", "less"))
  ylab <- list(eval)
  if(evalf.nme[i]=="meanLenOlapA"){
    df$obs <- log10(df$obs)
    evalf.nme[i] <- paste0("log10 ", evalf.nme[i])
  }
  if(evalf.nme[i]=="comOlap"){
    df$obs <- 100*df$obs/PERMTSUM$totlen[foi,"bin"]     
    evalf.nme[i] <- paste0(evalf.nme[i], "divTotBinlength")
  }
 
  # Order foi based on -log10(p-value)
  df$foi <- factor(df$foi, levels=c(df$foi[order(df$pval, decreasing=TRUE)]))
  
  plot.id <- paste0(out.id, "_", out.name, "_", evalf.nme[i])
  p <- ggplot(data=df, aes(y=foi, x=-log10(pval))) +
    geom_vline(xintercept=-log10(0.05), linetype="dashed", size=2) +
    geom_bar(stat="identity", aes(fill=alt)) + 
    geom_text(aes(label=format(obs, digits=3)), size=3, 
              position=position_stack(vjust=0.5)) + 
    labs(title=paste0(plot.id, "\n(", eval, ")"), y=NULL, 
         x=bquote(bold( "-log"["10"]~"(p-value)" ))
    ) +
    scale_fill_manual(values=col.v[c("greater", "less")%in%unique(df$alt)]) +
    bgr4 + 
    theme(axis.text.y=element_text(size=5, colour="black", face="bold",
                                   hjust=1),
          plot.title=element_text(size=5),
          panel.grid.major.y=element_line(colour="gray")) 
  
  ggsave(filename=paste0(out.dir, "/", plot.id, "_permplot.pdf"), 
         units="in", height=50, width=50, plot=p, limitsize=FALSE)
  print(paste0(eval, " plot done!"), quote=FALSE)
  rm(df, eval, p); gc()
  
}
# rm(list=ls())



