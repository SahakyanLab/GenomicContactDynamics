################################################################################
# Correlate C||align vs. C||kmer 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
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
csv.dir = paste0(wk.dir, "/out_compare_metric/label")
out.dir = paste0(wk.dir, "/out_compare_metric_pearsonplot")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
x <- list()
ijset.v = c("All", "LR", "SR") 
for(ijset in ijset.v){
  
  x[[ijset]] <- read.csv(file=paste0(csv.dir, "/", gcb, "_", ijset, "_pearsonR_alignVskmer.csv"),
                         header=TRUE)
  x[[ijset]][["chr"]] <- gsub(x=x[[ijset]][["chr"]], pattern="chr", replacement="")
  x[[ijset]] <- x[[ijset]][order(x[[ijset]][["PearsonR"]], decreasing=TRUE),]
  if(ijset=="All"){ chr.lev <- as.character(x[[ijset]][["chr"]]) }
  x[[ijset]] <- cbind(x[[ijset]], type=ijset)
  
}
x <- do.call("rbind", x)
x$chr <- factor(x$chr, levels=chr.lev); rm(chr.lev)
x$type <- factor(x$type, levels=c("All", "LR", "SR"))
rownames(x) <- NULL

p <- ggplot(data=x, aes(x=chr, y=PearsonR)) +
  geom_point(aes(col=type, shape=type, size=type)) +
  bgr2 +
  scale_shape_manual(values=c(8, 16, 16)) + 
  scale_size_manual(values=c(3, 2.5, 2.5)) + 
  theme(axis.text.x=element_text(size=18), 
        panel.grid.major.x=element_line(colour="gray80", linetype="dashed"))
ggsave(filename=paste0(out.dir, "/", gcb, "_pearsonR_alignVskmer_plot.pdf"),
       width=10, height=10)

# rm(list=ls()); gc()
