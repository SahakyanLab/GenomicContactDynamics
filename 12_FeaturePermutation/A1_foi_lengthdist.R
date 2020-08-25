################################################################################
# Length distribution of the features
# Mac, R/3.5.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
### OTHER SETTINGS #############################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib" 
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/reduced_b1b2b3b4")
# List of filenames of features of interest (refer to foi.dir)
# If foifile = NULL, all files in foi.dir
foifile = NULL #paste0(wk.dir, "/foifile/foifile_TF")
out.dir = paste0(wk.dir, "/out_foi_lengthdist")
### OTHER SETTINGS #############################################################
out.name = "foi_lengthdist_b1b2b3b4"
bin.len = 40000
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/finaliseFOI.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# List of features
foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)
foi.v.len <- length(foi.v)

toExport <- c("out.dir", "foi.v", "foi.dir")

p.lst <- list()
min.lst <- list()
for(foi in foi.v){

  # Only get the start and stop coordinates
  foi.bed <- read.table(file=paste0(foi.dir, "/", foi), stringsAsFactors=FALSE, 
                        header=FALSE)[,1:3]
  
  foi <- gsub(x=foi, pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
  df <- data.frame( length=log10(foi.bed[,3]-foi.bed[,2]+1L) )
  min.lst[[foi]] <- min(df$length)
  
  v <- round( c(min(df[,1]), max(df[,1])), digits=2)
  
  # Density plot
  p.lst[[foi]] <- ggplot(data=df, aes(x=length)) +
    geom_density(position="identity", fill="deepskyblue3", 
                 colour="deepskyblue3") + 
    scale_x_continuous(limits=c(1,10)) + 
    #scale_y_continuous(limits=c(0,80)) +
    geom_vline(xintercept=log10(bin.len), lty=2) + 
    labs(title=foi, 
         y=expression(bold( "Density" )), 
         x=bquote( bold("log"["10"]~"(L"^"foi"~")") ) 
    ) +
    bgr2 + 
    theme(plot.title=element_text(size=10))
  if(length(unique(v))!=1){
    p.lst[[foi]] <- p.lst[[foi]] + 
      annotate(geom="text", x=v, y=-0.1, size=7, label=v)
  }  
  
  rm(foi.bed); gc()
  
  print(paste0(foi, " done!"), quote=FALSE)
  
}

min.lst <- sort(unlist(min.lst, use.names=TRUE), decreasing=FALSE)
ind <- match(names(min.lst), names(p.lst))

p.lst <- ggarrange(plotlist=p.lst[ind], nrow=3, ncol=6,
                   common.legend=TRUE, legend="right")
ggexport(p.lst, width=52, height=24,
         filename=paste0(out.dir, "/", out.name, ".pdf"))

# rm(list=ls())
