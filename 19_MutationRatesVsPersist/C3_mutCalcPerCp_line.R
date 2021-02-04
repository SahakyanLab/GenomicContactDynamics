################################################################################
# Plot median and mean of mutation calculations per Cp. Calculate significant
# relative to Cp=1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    binmx.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/binmx/out_bindata_1perc_HiCNorm"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
mutCp.dir = paste0(wk.dir, "/out_mutCalcPerCp/WT_SEQ_rowSum")
out.dir = paste0(wk.dir, "/out_mutCalcPerCp_line/WT_SEQ_rowSum")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
src.id = "hg38ToHg19" # "Hg19" | "hg38ToHg19"
mut.v = c("All", "A>C", "A>G", "A>T", "C>A", "C>G", "C>T")
calc.v = c("Tmut", "Nmsite", "TmutDIVNmsite", "Nmsitenorm")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
p.lst <- NULL
for(mut in mut.v){
  
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=TRUE)
  load(file=paste0(mutCp.dir, "/", mut.id, "_", src.id, "_mutCalcPerCp.RData"))
 
  #-------------------Boxplots
  x <- reshape2::melt(data=MUTCP.DF[,-(1:2)], id.vars=c("mutbin", calc.v))
  rm(MUTCP.DF); gc()
  # Take only bins with Cp >= 1
  x <- x[x$value==1,colnames(x)!="value"]
  colnames(x)[colnames(x)=="variable"] <- "Cp"
  Cp.v <- sort(as.numeric(unique(x$Cp)))
  
  for(calc in calc.v){
    
    mean.df <- aggregate(x=x[[calc]], by=list(x$Cp), FUN=mean, na.rm=TRUE)
    mean.df$Group.1 <- as.numeric(mean.df$Group.1)
    mean.df$id <- "MEAN"
    rownames(mean.df) <- paste0(mean.df$Group.1, "_MEAN")
    med.df <- aggregate(x=x[[calc]], by=list(x$Cp), FUN=median, na.rm=TRUE)
    med.df$Group.1 <- as.numeric(med.df$Group.1)
    med.df$id <- "MEDIAN"
    rownames(med.df) <- paste0(med.df$Group.1, "_MEDIAN")
    df <- rbind.data.frame(mean.df, med.df)
    
    colnames(df) <- c("Cp", "value", "id")
    rm(mean.df, med.df)
    
    # Wilcoxon-Mann-Whitney test (non-parametric, unpaired/independent groups) 
    # to compare orig and shuff values per Cp
    # The test assumes that the shape of the two distributions are similar so
    # check the boxplots
    # If both x and y are given and paired is FALSE, a Wilcoxon rank sum test
    # (equivalent to the Mann-Whitney test) is carried out.
    df$pval <- df$alt <- NA
    
    for(Cp in Cp.v){
      
      # pos --> Cp==N > Cp==1
      # 0   --> Cp==N = Cp==1
      # neg --> Cp==N < Cp==1
      id.mean <- paste0(as.character(Cp), "_MEAN")
      id.med <- paste0(as.character(Cp), "_MEDIAN")
      
      alt.v <- c(MEAN=df[id.mean,"value"]-df["1_MEAN","value"],
                 MEDIAN=df[id.med,"value"]-df["1_MEDIAN","value"])
      alt.v <- alt.v/abs(alt.v)
      if( !all(alt.v%in%c(-1, 1, NaN)) ){
        stop("Checkpoint 1.")
      }
      
      alt.v <- c(`-1`="less", `1`="greater","NaN"="two.sided")[as.character(alt.v)]
      names(alt.v) <- c("MEAN", "MEDIAN")
      
      # The one-sided alternative "greater" is that x is shifted to the right of y
      df[id.mean,"pval"] <- wilcox.test(x=x[x$Cp==Cp,calc], y=x[x$Cp==1,calc],
                                               alternative=alt.v["MEAN"], paired=FALSE)$p.value 
      df[id.med,"pval"] <- wilcox.test(x=x[x$Cp==Cp,calc], y=x[x$Cp==1,calc],
                                       alternative=alt.v["MEDIAN"], paired=FALSE)$p.value 
      
      df[id.mean,"alt"] <- alt.v["MEAN"] 
      df[id.med,"alt"] <- alt.v["MEDIAN"] 
      
      rm(alt.v)
      
    } # Cp.v for loop end
    
    # Plot
    df$pval <- c(`0`="ns", `1`="sig")[ as.character(as.numeric(df$pval<0.05)) ]
    df$pval <- factor(df$pval, levels=c("ns", "sig"))
    df$Cp <- factor(df$Cp, levels=c(as.character(sort(unique(df$Cp)))))
    
    id <- paste0(gcb, "_", src.id, "_", mut, "_", calc)
    p.lst[[id]] <- ggplot(data=df, aes(x=Cp, y=value, group=id)) +
      geom_point(aes(colour=id, alpha=pval, shape=alt), size=3) +
      scale_alpha_manual(values=c(0.3,1)) + 
      labs(x="Cp", y=calc, title=id) + 
      bgr2

  }  # calc.v for loop end
  
  print(paste0(mut, " done!"), quote=FALSE)
  
} # mut.v for loop end

p.arr <- ggarrange(plotlist=p.lst, nrow=length(mut.v), ncol=length(calc.v), legend=NULL)
ggexport(p.arr, height=length(mut.v)*10, width=length(calc.v)*10,
         filename=paste0(out.dir, "/", gcb, "_", src.id, "_scatplot.pdf"))

# rm(list=ls()); gc()

