################################################################################
# Arcplot distinguishing between long- and short-range persistent contacts.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/17_Structure"
    data.dir= "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/17_Structure"
    data.dir= "/t1-data/user/ltamon/Database"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/17_Structure"
    data.dir= "/home/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_arcplotB")
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
# Expands warnings
#options(warn=1)
gcb = "min2Mb"
chr.v = paste0("chr", c(22:1, "X"))
bin.len = 40000
# topCp = 3 would mean the top 3 Cps so 19:21
topCP = 3
# Range of Cp scores
cp.v = 21:1
gap.perc.thresh = 20
# xlim max value
xlim = NULL
plotOnly = FALSE
showSRij = TRUE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(R4RNA)
library(RColorBrewer)
library(grDevices)
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Chromosome length file
chrLen.df <- read.table(file=chrLenfile, as.is=FALSE, header=TRUE,
                        colClasses=c("character", "integer", "integer"))
cp <- cp.v[1:topCP]
cp <- sort(cp, decreasing=FALSE)
for(chr in chr.v){
  out.name <- paste0(chr, "_topCP", topCP, "_gapPerc", gap.perc.thresh, "_arcLSrange")
  if(plotOnly==FALSE){
    df <- list()
    # Load PERSIST.MX
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    incl.TF <- PERSIST.MX$ntis%in%cp
    if(sum(incl.TF)==0){
      print(paste0(chr, " skipped!"), quote=FALSE)
      next
    }
    gap.perc <- ((PERSIST.MX$hits$j[incl.TF]-PERSIST.MX$hits$i[incl.TF])*bin.len)-bin.len
    gap.perc <- 100*gap.perc/( chrLen.df[chrLen.df$chromosome==chr,"length.bp"] )
    df$bottom <- cbind(PERSIST.MX$hits[incl.TF, c("i", "j")], length=as.integer(1), 
                       value=PERSIST.MX$ntis[incl.TF], gap.perc=gap.perc)
    rm(PERSIST.MX, gap.perc, incl.TF);gc()
    #---------------------------------------
    # Select long-range persistent contacts based on gap threshold
    top.TF <- df$bottom$gap.perc>gap.perc.thresh
    if(sum(top.TF)==0){
      print(paste0(chr, ":No long-range contacts (top panel)"), quote=FALSE)
      df$top <- NULL
    } else {
      # Long-range persistent contacts
      df$top <- df$bottom[top.TF,]
      df$top <- colourByValue(as.helix(df$top), breaks=c(0,cp), right=TRUE,
                              include.lowest=FALSE, get=TRUE, 
                              cols=adjustcolor(tail(brewer.pal(n=4, name="Reds"), n=topCP), alpha.f=0.8))
      # Short-range persistent contacts
      df$bottom <- df$bottom[!top.TF,]
    }
    df$bottom$col <- adjustcolor("gray50", alpha.f=0.4)
    save(df, file=paste0(out.dir, "/", out.name, ".RData"))
  } else {
    file.nme <- paste0(out.dir, "/", out.name, ".RData")
    if( file.exists(file.nme) ) {
      load(file=file.nme)
    } else {
      print(paste0(chr, " skipped!"), quote=FALSE)
      next
    }
  }
  #---------------------------------------
  if( !is.null(xlim) ){
    v <- df$bottom[1,]
    v[1] <- v[2] <- xlim
    v[3] <- 0
    if( !is.null(df$top) ){
      df$top <- rbind(df$top, v)
    }
    df$bottom <- rbind(df$bottom, v); rm(v)
  }
  if(showSRij==FALSE){
    df$bottom <- NULL
    out.name <- paste0(out.name, "_LRonly")
  }
  #---------------------------------------
  pdf(file=paste0(out.dir, "/", out.name, ".pdf"), 
      width=10, height=10)
  plotDoubleHelix(as.helix(df$top), as.helix(df$bottom), shape="circle", line=TRUE, 
                  arrow=FALSE, lwd=2, oma=c(0,0,0,0))
  legend("bottomright", legend=sort(cp, decreasing=FALSE), 
         fill=attr(df$top, "fill"), inset=0.1, bty="n", border=NA, 
         cex=1, title=expression(bold("c"["p"]))
         )
  top.ij <- sum(df$top$length==1)
  bot.ij <- sum(df$bottom$length==1)
  mtext(paste0(chr, "_", gcb, "_topCP=", topCP,  "_gap>", gap.perc.thresh, 
               "%_totij=", top.ij+bot.ij, "_topij=", top.ij, "_botij=", bot.ij), 
        side=3, line=-1, adj=0.5, cex=0.5)
  dev.off()
  
  rm(df); gc()
  print(paste0(chr, " done!"), quote=FALSE)
}

# rm(list=ls())