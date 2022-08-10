################################################################################
# Plot number of contacts formed by each bin per chr
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(wk.dir, "/out_basePersist")
out.dir = paste0(wk.dir, "/out_binVsNij")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
topCP = 3; cp.v = 1:21
# If ct = "All", no filtering based on cell/tissue.
ct = "All"; ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
                     "AG", "Ov", "Bl", "MesC","MSC", "NPC", "TLC", "ESC", "FC", "LC")
# Gap in terms of bin
gap.val.v = c(50, 75, 100, 125, 150, 160, 180, 200, 250, 300, 350, 400, 600, 
              650, 700)
lab.top = 20
Nij.thresh = 3
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(ggrepel)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
cp.v <- sort(unique(cp.v), decreasing=FALSE)
ct.id <- NULL
if( !is.null(ct) ){ ct.id <- ct; ct.id <- paste0(ct.id[!is.na(ct.id)], "_") }

mx <- matrix(data=NA, nrow=length(chr.v), ncol=length(gap.val.v), 
             dimnames=list(chr.v, gap.val.v))

for(chr in chr.v){
  
  load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, "_topCP4.RData"))
  max.bin <- max(unique(c(PERSIST.MX$hits$i, PERSIST.MX$hits$j)))
  incl.TF <- PERSIST.MX$ntis%in%tail(x=cp.v, n=topCP)
  if(ct%in%ct.v){
    incl.TF <- incl.TF & PERSIST.MX$hits[[ct]]>0
    print("Filtering based on cell/tissue.", quote=FALSE)
  } else if(ct=="All"){
    print("No filtering based on cell/tissue.", quote=FALSE)
  } else {
    stop("Invalid ct argument.")
  }
  PERSIST.MX$hits <- PERSIST.MX$hits[incl.TF,c("i", "j")]; rm(incl.TF)
  PERSIST.MX$ntis <- NULL
  gap.v <- PERSIST.MX$hits$j-PERSIST.MX$hits$i
  
  p.lst <- list()
  out.name <- paste0(gcb, "_", chr, "_topCP", topCP, "_", ct)
  for(gap.val in gap.val.v){
    
    bins.df <- table(unlist(PERSIST.MX$hits[gap.v>gap.val,c("i","j")], use.names=FALSE))
    bins.df <- as.data.frame(cbind(Bin=as.numeric(names(bins.df)), Nij=unname(bins.df)))
    bins.df$lab <- bins.df$Nij
    if(nrow(bins.df)>0){
      bins.df$lab <- NA
      bins.df <- bins.df[order(-bins.df$Nij),]
      lab.lim <- min(nrow(bins.df),lab.top)
      bins.df$lab[1:lab.lim] <- bins.df$Bin[1:lab.lim]
      thresh.TF <- bins.df$Nij>Nij.thresh
      if(sum(thresh.TF)>0){
        mx[chr,as.character(gap.val)] <- paste(bins.df[thresh.TF,"Bin"], collapse=";")
      }
    }
    
    p.lst[[paste0(chr, "_", as.character(gap.val))]] <- ggplot(data=bins.df, 
                                                               aes(x=Bin, y=Nij, label=lab)) +
      geom_text_repel(size=2, segment.colour=adjustcolor("gray",0.5), 
                      na.rm=TRUE) + 
      geom_point(size=2, shape=1, colour="red") + 
      scale_y_continuous( limits=c(0,150), breaks=sort(c(3, 5, seq(0, 150, by=10))) ) + 
      scale_x_continuous(limits=c(1,max.bin)) + 
      labs(title=paste0(out.name, "_gapBin", gap.val, "_totNij=", sum(bins.df$Nij)/2)) + 
      bgr2 + 
      theme(panel.grid.major.y=element_line(colour=adjustcolor("blue",0.5)))
    
    rm(bins.df); gc()
  } 
  rm(PERSIST.MX, max.bin, gap.v); gc()
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=3, ncol=5)
  ggexport(p.arr, height=30, width=50,
           filename=paste0(out.dir, "/", out.name, ".pdf"))
  rm(p.lst, p.arr); gc()
  print(paste0(chr, " done!"), quote=FALSE)
  
}
write.csv(cbind(chr=rownames(mx),mx), 
          file=paste0(out.dir, "/", gcb, "_topCP", topCP, "_", ct, "_Nijgrthan", Nij.thresh, ".csv"),
          row.names=FALSE, quote=FALSE)

# rm(list=ls()); gc()



