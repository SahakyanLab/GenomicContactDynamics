################################################################################
# Correlate C||align vs. C||kmer 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/11_Complementarity"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir = "/t1-data/user/ltamon/Database"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir = "/home/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
compl.dir = paste0(wk.dir, "/out_constraints/merged_final")
#compl.dir = paste0(wk.dir, "/out_constraints")
out.dir = paste0(wk.dir, "/out_compare_metric/label")
### OTHER SETTINGS #############################################################
chr.v = paste("chr", c(1:22, "X"), sep="")
ijset = "All" # Contact set, All | LR (long-range) | SR(short-range)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
#source(paste0(lib, "/lmEqn_string.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
gcb = "min2Mb"
p.lst <- list()
R.v <- rep(NA, times=length(chr.v))
names(R.v) <- chr.v
for(chr in chr.v){
  
  df <- NULL
  # Load CII.MX kmer
  load(file=paste0(compl.dir, "/", chr, "_kmer_", gcb, ".RData"))
  
  # Contact set, All | long-range | short-range
  if(ijset!="All"){
    
    if(ijset=="LR"){
      
      ij.TF <- !is.na(CII.MX[,"Cp"])
      print("Long-range contacts only...", quote=FALSE)
      
    } else if(ijset=="SR"){
      
      ij.TF <- is.na(CII.MX[,"Cp"])
      print("Short-range contacts only...", quote=FALSE)
      
    } else {
      stop("Wrong ijset argument!")
    }
    
  } else {
    
    ij.TF <- rep(TRUE, times=nrow(CII.MX))
    print("All contacts...", quote=FALSE)
    
  } 
  
  df$kmer <- CII.MX[,"C||"]
  # Total contacts, same for align
  len0 <- nrow(CII.MX)
  len <- sum(ij.TF)
  # Contacts with kmer C||
  nonNA.TF <- !is.na(df$kmer)
  rm(CII.MX); gc()
  
  # Load CII.MX align
  load(file=paste0(compl.dir, "/", chr, "_align_", gcb, ".RData"))
  if(len0!=nrow(CII.MX)){ stop("kmer and align CII.MX have different dim.") }
  df$align <- CII.MX[,"C||"]
  
  # Contacts with kmer and align C|| and the right contact type
  nonNA.TF <- nonNA.TF & !is.na(df$align) 
  rm(CII.MX, len0); gc()
  
  df <- do.call("cbind.data.frame", df)
  # df contains values for contacts with both types of scores (no NAs)
  df <- df[nonNA.TF & ij.TF,]
  rm(nonNA.TF, ij.TF); gc()
  
  tot.ij <- format(as.numeric(len), scientific=TRUE, digits=4)
  nonNA.ij <- format( nrow(df)/len*100, scientific=TRUE, digits=4)
  
  R.v[chr] <- cor(x=df[,"kmer"], y=df[,"align"], method="pearson")
  
  p.lst[[chr]] <- ggplot(data=df, aes(x=kmer, y=align)) +
    geom_point(fill=adjustcolor("gray59", alpha.f=0.2), color=adjustcolor("gray59", alpha.f=0.2)) +
    stat_smooth(geom="smooth", method="lm", formula=y~x, se=TRUE, 
                n=80, span=0.75, level=0.95, aes(colour="darkred"), size=3) +
    guides(colour=FALSE) + 
    labs(
      title=paste0(chr, "_", gcb, "_", tot.ij, ijset, "ij_", nonNA.ij, "nonNACII_"), 
      x=bquote(bold( "C"["||"]~"kmer" )), 
      y=bquote(bold( "C"["||"]~"align" )) 
    ) + 
    #labs(title=NULL, x=NULL, y=NULL) +
    #theme(axis.text.x=element_blank(), axis.text.y=element_blank()) +
    bgr3 #+
    #annotate(geom="text", x=min(df$kmer), y=c(max(df$align), max(df$align)-0.02),
    #annotate(geom="text", x=min(df$kmer), y=max(df$align),
    #         size=8, parse=TRUE, hjust=0, vjust=1,
    #         label=lmEqn_string(x="kmer", y="align", data=df)
    #) 
  
  ggsave(filename=paste0(out.dir, "/", chr, "_", gcb, "_", ijset, "_metricCor.png"),
         width=10, height=10, units="in", dpi=1200, plot=p.lst[[chr]])
  
  rm(df); gc()
  
}

R.v <- stack(R.v)
colnames(R.v) <- c("PearsonR", "chr")
R.v$chr <- factor(as.character(R.v$chr), levels=as.character(R.v$chr))
write.csv(R.v, file=paste0(out.dir, "/", gcb, "_", ijset, "_pearsonR_alignVskmer.csv"), row.names=FALSE)

p <- ggplot(data=R.v, aes(x=chr, y=PearsonR)) +
  #geom_line(size=2.5) +
  geom_point(size=5, colour="darkred") +
  labs(title=paste0(gcb, "_pearsonR_alignVskmer_", ijset), x="chr", 
       y=expression(paste("Pearson's")~italic("r"))) +
  bgr2 
  
ggsave(filename=paste0(out.dir, "/", gcb, "_", ijset, "_alignVskmer_PCoeff.pdf"),
       width=10, height=10, plot=p)

p.arr <- ggarrange(plotlist=p.lst, nrow=4, ncol=6,
                   legend=NULL)
#ggexport(p.arr, width=6000, height=4000, units="px",
ggexport(p.arr, width=4500, height=3000, units="px",
#ggexport(p.arr, width=3000, height=2000, units="px",
         filename=paste0(out.dir, "/", gcb, "_", ijset, "_metricCor.png"))

# rm(list=ls()); gc()
