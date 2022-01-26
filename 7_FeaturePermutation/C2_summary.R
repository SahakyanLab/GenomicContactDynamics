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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
Cp.id = "CpAllCs1perc" #"Cp21" #"CptopCP3" #CpAllCs1perc
permt.dir = paste0(wk.dir,"/z_ignore_git_COMPLETE/out_associate_", Cp.id)

out.dir = paste0(wk.dir,"/out_summary/feat_844_raw")
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/raw_ALL_associated")

#out.dir = paste0(wk.dir,"/out_summary/feat_526_raw")
#foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/raw_associated")

foifile = NULL #paste0(wk.dir, "/foifile/foifile_test1")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
permt.id = "nperm10000_seed662_mxmskfr0"
eval.f.v = c("numOlapA", "numOlapAwithin", "numOlapB", "numOlapBwithin",
             "comOlap", "meandist", "meanLenOlapB")
A = "bin"
B = "foi"
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
evalf.v.len <- length(eval.f.v)
out.id <- paste0(permt.id, "_", Cp.id)

permt.v <- list.files(path=permt.dir, pattern=permt.id)
permt.v <- permt.v[grepl(x=permt.v, pattern=paste0(Cp.id, "_permtest.RData"), 
                         fixed=TRUE)]

if(plotOnly==FALSE){
 
  foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)
  ct.v <- unlist(
    lapply(X=strsplit(x=foi.v, split="ct_|\\_"), FUN=function(x)x[2])
  )
  names(ct.v) <- foi.v
  foilab.v <- sapply(X=foi.v, simplify=TRUE, USE.NAMES=FALSE, FUN=function(foi){
    ct <- ct.v[foi]
    pat <- paste(paste0("ct_", ct, "_"), "foi\\_|desc\\_|\\.bed", sep="|")
    foi <- gsub(x=foi, pattern=pat, replacement="")
    foi <- paste0(ct, "_", foi)
    return(foi)
  }); rm(ct.v)
  foi.v.len <- length(foi.v)
  
  # Initialize output tables
  PERMTSUM <- list()
  PERMTSUM$pval <- matrix(data=NA, ncol=length(eval.f.v), nrow=foi.v.len,  
                          dimnames=list(foilab.v, eval.f.v))
  PERMTSUM$obs <- PERMTSUM$pval
  PERMTSUM$alt <- PERMTSUM$pval
  PERMTSUM$numANDlen <- matrix(data=NA, ncol=4, nrow=foi.v.len, 
                               dimnames=list(foilab.v, c("Anum", "Alen", "Bnum", "Blen")))
  num.v  <- dimnames(PERMTSUM$numANDlen)[[2]]
  PERMTSUM$foifile <- foi.v; rm(foi.v)
                         
  for(i in 1:foi.v.len){
    
    foi <- foilab.v[i]
    fle <- permt.v[grepl(x=permt.v, pattern=paste0(gcb, "_", foi, "_"), fixed=TRUE)]
    if(length(fle)!=1){ stop("Multiple file.") }
    load(file=paste0(permt.dir, "/", fle))
    
    for(x in 1:evalf.v.len){
      eval <- eval.f.v[x]
      PERMTSUM$pval[foi,eval] <- as.numeric(PERMT[[eval]]$pval)
      PERMTSUM$obs[foi,eval]  <- as.numeric(PERMT[[eval]]$observed)
      PERMTSUM$alt[foi,eval]  <- PERMT[[eval]]$alternative 
      if(x==1){
        PERMTSUM$numANDlen[foi,num.v] <- PERMT$region[num.v]
      }
    }
    rm(PERMT); gc()
    print(paste0(foi, " done!"), quote=FALSE)
    
  }
  save(PERMTSUM, file=paste0(out.dir, "/", out.id, "_permtsum.RData"))
  
} else {
  load(file=paste0(out.dir, "/", out.id, "_permtsum.RData"))
}

eval.f.v <- colnames(PERMTSUM$pval)
for(i in 1:evalf.v.len){
  
  evalf.nme <- eval <- eval.f.v[i]
  df <- data.frame(foi=rownames(PERMTSUM$pval), obs=PERMTSUM$obs[,eval], alt=PERMTSUM$alt[,eval],
                   pval=PERMTSUM$pval[,eval], row.names=NULL, stringsAsFactors=FALSE)
 
  df$alt <- factor(df$alt, levels=c("greater", "less"))
  ylab <- list(eval)
  if( eval%in%c("meanLenOlapA", "meanLenOlapB") ){
    df$obs <- log10(df$obs)
    evalf.nme <- paste0("log10", evalf.nme)
  }
  if(eval=="comOlap"){
    bin.letter <- c("A","B")[c(A,B)=="bin"]
    eval(parse(text=paste0(
      "df$obs <- 100*df$obs/PERMTSUM$numANDlen[,'", bin.letter, "len']"   
    )))
    evalf.nme <- paste0(evalf.nme, "divTotBinlength")
  }
  if( eval%in%c("numOlapA", "numOlapAwithin") ){
    df$obs <- 100*df$obs/PERMTSUM$numANDlen[,"Anum"]
    evalf.nme <- paste0(evalf.nme, "divTotAnum")
  }
  if( eval%in%c("numOlapB", "numOlapBwithin") ){
    df$obs <- 100*df$obs/PERMTSUM$numANDlen[,"Bnum"]
    evalf.nme <- paste0(evalf.nme, "divTotBnum")
  }
  if( eval%in%c("numOlapA", "numOlapAwithin", "numOlapB", "numOlapBwithin", "comOlap") ){
    if( any(df$obs>100) ){
      stop(paste0(eval, ": Percentage > 100%"))
    }
  }
  
  # Order foi based on -log10(p-value)
  df$foi <- factor(df$foi, levels=c(df$foi[order(df$pval, decreasing=TRUE)]))
  
  plot.id <- paste0(out.id, "_", evalf.nme)
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

# rm(list=ls()); gc()



