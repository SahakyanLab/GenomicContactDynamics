################################################################################
# Plot homer results
# deva, R/3.5.0-newgcc, gcc/4.9.2
# Mac, R/3.5.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/9_Motif"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/15_Motif"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
results.dir = paste0(wk.dir, "/z_ignore_git/out_motif/chrALL_uniquebins_targCp21_bgrCp1_seed845_sampsize0.7")
### OTHER SETTINGS #############################################################
topN.v <- c(50,50) # c(<known>, <deva>)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(marge)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
res.v <- c("known", "denovo")
names(topN.v) <- res.v
exists.v <- file.exists(paste0(results.dir, "/knownResults.txt"), paste0(results.dir, "/homerMotifs.all.motifs"))
size.v <- file.info(paste0(results.dir, "/knownResults.txt"), paste0(results.dir, "/homerMotifs.all.motifs"))$size
size.v <- size.v>0
res.v <- res.v[exists.v & size.v]
stopifnot(length(res.v)>=0)

for( res in res.v ){
  
  topN <- topN.v[res]
  
  if(res=="known"){
    
    tbl <- fread(file=paste0(results.dir, "/knownResults.txt"), data.table=FALSE)
    topN <- ifelse(topN > nrow(tbl), nrow(tbl), topN)
    tbl <- tbl[order(tbl$`Log P-value`, decreasing=FALSE),][1:topN,]
    # Format table
    targPerc <- as.numeric(unlist(
      strsplit(x=tbl$`% of Target Sequences with Motif`, split="%")
    ))
    bgrPerc <- as.numeric(unlist(
      strsplit(x=tbl$`% of Background Sequences with Motif`, split="%")
    )) 
    name <- strsplit(x=tbl$`Motif Name`, split="\\/")
    name <- paste(unlist(lapply(X=name, FUN=function(x) x[1]), use.names=FALSE),
                  tbl$`% of Background Sequences with Motif`, sep="_")
    neglog10pval <- -(log10(tbl$`P-value`))
    
    # Get number of target and bgr sequences used and add to title
    val <- colnames(tbl)[c(6,8)]
    val <- strsplit(x=val, split="\\(of |\\)")
    val <- unlist(lapply(X=val, FUN=function(x) x[2]))
    plottitle <- paste0(rev(strsplit(x=results.dir, split="\\/")[[1]])[1],
                        "_", val[1], "targSeq_", val[2], "bgrSeq")
    
  } else if(res=="denovo"){
    
    tbl <- data.frame(read_denovo_results(path=results.dir, homer_dir=TRUE))
    topN <- ifelse(topN > nrow(tbl), nrow(tbl), topN)
    tbl <- tbl[order(tbl$log_p_value_detection, decreasing=FALSE),][1:topN,]
    # Format table
    targPerc <- as.numeric(tbl$tgt_pct)*100
    bgrPerc <- as.numeric(tbl$bgd_pct)*100
    name <- paste(tbl$consensus, paste(bgrPerc, "%", sep=""), sep="_")
    neglog10pval <- -(log10(exp(tbl$log_p_value_detection)))
    
  } else {
    stop("Checkpoint.")
  }
  
  # Generic
  name <- paste(name, 1:length(name), sep="-")
  fc <- log2(targPerc/bgrPerc)
  
  # Make df for plotting
  df <- data.frame(name=name, FC=fc, neglog10pval=neglog10pval)
  # Plot by increasing p-value
  df$name <- factor(x=as.factor(df$name), levels=rev(as.character(df$name)))
  
  rm(tbl, name, targPerc, bgrPerc, neglog10pval, topN)
  
  p <- ggplot(data=df, aes(x=neglog10pval, y=name)) +
    geom_point(aes(colour=FC), size=4) +
    geom_vline(xintercept=-log10(0.05), linetype="dashed", colour="red", size=1) +
    scale_colour_gradient2(high="red4", mid="mistyrose", low="blue", midpoint=0, 
                           space="Lab", na.value="gray50", guide="colourbar", 
                           aesthetics="colour") + 
    labs(title=paste0(plottitle, "_", res),
         x="-log10(p-value)", y=NULL, colour="FC targ/bgr") + 
    theme(panel.grid.major.y=element_line(colour="gray"),
          panel.grid.minor.y=element_line(colour="gray"),
          panel.grid.major.x=element_line(colour="gray"),
          panel.background=element_rect(colour="gray22", 
                                        size=1, fill="ghostwhite"),
          plot.title=element_text(size=8, colour="black"),
          axis.title.y=element_blank(),
          axis.text.y=element_text(size=10, colour="black", face="bold"),
          axis.title.x=element_text(size=10, colour="black", face="bold"),
          axis.text.x=element_text(size=10, colour="black", face="bold"),
          legend.text=element_text(size=10, face="bold"),
          legend.title=element_text(size=10, face="bold")
    )
  
  ggsave(filename=paste0(results.dir, "/result_plot_", res, ".pdf"), plot=p,
         width=9, height=10)
  
  rm(p); gc()
  
  print(paste0(res, " done!"), quote=FALSE)
  
}

# rm(list=ls())
