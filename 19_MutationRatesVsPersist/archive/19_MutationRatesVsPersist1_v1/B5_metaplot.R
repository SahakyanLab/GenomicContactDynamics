################################################################################
# Two plots:
# a. Get mean %bin for each 40kb-bin then calculate fold change relative to mean of 
# actual HiC bin. Plot is log2(FC) X midpoint of each bin(x10^4)
# b. See what feature go down or up with increasing Cp. Per feature, plot
# mean %bin of actual HiC bin Vs contact persistence. 
# deva, R/3.6.0-newgcc, gcc/4.9.2

# ERROR: missing points in plot, this is due to preset y-axis min and max limits 
# in plotParam (i.e. -1, 1), change if necessary
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    data.dir = "/Users/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
foifile = paste0(wk.dir, "/foifile/foifile_mut")
# File of feature grouping
featgrpfile = paste0(wk.dir, "/features_group")
fetacp.dir = paste0(wk.dir, "/out_combChr_count")
out.dir = paste0(wk.dir, "/out_metaplot")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
out.id = "mutation"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(reshape)
library(ggplot2)
library(ggsci)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
myplot <- function( df=df, title=paste0(foilab, "_Cp=", Cp) ){
  
  p <- ggplot(data=df, aes(x=bin, y=value)) +
    geom_col(col=col.v[foi.v[f]]) + 
    labs(title=title, x="Position", y="# Mut") + 
    theme( axis.text.x=element_text(size=10) ) +
    bgr2
  return(p)
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
foi.v <- readLines(con=foifile)
foi.v.len <- length(foi.v)

col.v <- ggsci::pal_npg(palette="nrc", alpha=0.5)(foi.v.len)
names(col.v) <- foi.v

for(f in 1:foi.v.len){
  
  # Load FETACP.MX
  load(file=paste0(fetacp.dir, "/chrALL_", gcb, "_", 
                   gsub(x=foi.v[f], pattern="ct\\_|foi\\_|desc\\_|\\.bed", 
                        replacement=""), "_fetacp.RData"))
  Cp.v <- as.numeric(rownames(FETACP.MX))
  
  foilab <- gsub(x=foi.v[f], pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
  
  if(f==1){
    pos.v <- as.numeric(colnames(FETACP.MX))
  }
  
  df <- melt(FETACP.MX)
  rm(FETACP.MX); gc()
  colnames(df) <- c("Cp", "bin", "value")
  df$Cp <- as.numeric(df$Cp)
  
  p <- myplot(df=aggregate(formula=value~bin, data=df, FUN=sum), 
              title=foilab)
  ggsave(filename=paste0(out.dir, "/chrALL_", gcb, "_", foilab, ".pdf"), 
         width=10, height=10, plot=p)
  
  p.lst <- list()
  for(Cp in Cp.v){
    
    temp <- aggregate(formula=value~bin, data=df[df$Cp==Cp,], FUN=sum)
    id <- paste0(foilab, "_Cp=", Cp)
    p.lst[[id]] <- myplot(df=temp, title=id)
    rm(temp)

  }

  p.arr <- ggarrange(plotlist=p.lst, nrow=3, ncol=7, legend=NULL)
  ggexport(p.arr, height=30, width=70, 
           filename=paste0(out.dir, "/chrALL_", gcb, "_", foilab, "_perCp.pdf"))
  
  rm(p, foilab, df, p.arr, p.lst); gc()
  
} # foi.v for loop end

# rm(list=ls()); gc()



