################################################################################
# Arcplot distinguishing between long- and short-range persistent contacts.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/13_ArcPlot"
    data.dir= "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/17_Structure"
    data.dir= "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_arcplotB")
chrLenfile = paste0(wk.dir, "/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
# Expands warnings
options(warn=1)

gcb = "min2Mb"
bin.len = 40000
chr.v = paste0("chr", c(22:1, "X"))

# topCP=3 would means the top 3 values of cp.v; topCP=-3 means the bottom 3 values
# of cp.v
topCP = 3; cp.v = 1:21
ct = "FC"; ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
                    "AG", "Ov", "Bl", "MesC","MSC", "NPC", "TLC", "ESC", "FC", "LC")
# Gap in terms of % of chr length ("Perc") or bin ("Bin")
gap.type = "Bin" 
gap.val = 150

# xlim max value to make x-axis consistent across chr
xlim = NULL
plotOnly = FALSE
# If showSRij = TRUE, contacts below gap threshold will be shown at the bottom
# panel
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
#---------------------------------------
cp.v <- sort(cp.v, decreasing=FALSE)
if(topCP<0){ cp.v <- rev(cp.v) }
cp <- sort(tail(cp.v, n=abs(topCP)), decreasing=FALSE)
#---------------------------------------
ct.id <- NULL
if( !is.null(ct) ){ ct.id <- ct; ct.id <- paste0(ct.id[!is.na(ct.id)], "_") }
#---------------------------------------
for(chr in chr.v){
  
  out.name <- paste0(gcb, "_", chr, "_", ct.id, "topCP", topCP, "_gap", gap.type, 
                     gap.val, "_arcLSrange")
  
  if(plotOnly==FALSE){
    
    df <- list()
    
    # Load PERSIST.MX
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    if( is.null(ct) ){
      incl.TF <- PERSIST.MX$ntis%in%cp 
    } else if(ct%in%ct.v){
      incl.TF <- PERSIST.MX$ntis%in%cp & PERSIST.MX$hits[[ct]]>0
    } else {
      stop("Invalid ct argument. It should be NULL or an element of ct.v")
    }
    ij.df <- cbind(PERSIST.MX$hits[incl.TF,c("i", "j")], Cp=PERSIST.MX$ntis[incl.TF])
    rm(PERSIST.MX, incl.TF); gc()
    if(nrow(ij.df)==0){
      print(paste0(chr, " skipped!"), quote=FALSE)
      next
    }
    
    gap.bin.v <- (ij.df$j-ij.df$i)-1
    if(gap.type=="Perc"){
      chr.len <- chrLen.df[chrLen.df$chromosome==chr,"length.bp"]
      gap.TF <- (100*(gap.bin.v*bin.len/chr.len)) >= gap.val
    } else if(gap.type=="Bin"){
      gap.TF <- gap.bin.v >= gap.val
    } else {
      stop("Invalid gap.type. It should be Bin or Perc only.") 
    }
    rm(gap.bin.v)
    #---------------------------------------
    # Select long-range persistent contacts based on gap threshold
    if(sum(gap.TF)==0){
      print(paste0(chr, ":No long-range contacts passing the gap threshold (top panel)"), 
            quote=FALSE)
      df$top <- NULL
    } else {
      # Longest of the long-range persistent contacts
      df$top <- cbind(ij.df[gap.TF, c("i", "j")], length=1L, value=ij.df$Cp[gap.TF])
      df$top <- colourByValue(as.helix(df$top), breaks=c(0,cp), right=TRUE,
                              include.lowest=FALSE, get=TRUE, 
                              cols=adjustcolor(tail(brewer.pal(n=4, name="Reds"), 
                                                    n=topCP), alpha.f=0.8))
      # Rest of long-range persistent contacts
      bottom.TF <- !gap.TF
      df$bottom <- cbind(ij.df[bottom.TF, c("i", "j")], length=1L, value=ij.df$Cp[bottom.TF])
      rm(bottom.TF)
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
  legend( "bottomright", legend=sort(cp, decreasing=FALSE), 
          fill=attr(df$top, "fill"), inset=0.1, bty="n", border=NA, 
          cex=1, title=expression(bold("c"["p"])) )
  # I did it this way to not count extra row added by xlim
  top.ij <- sum(df$top$length==1)
  bot.ij <- sum(df$bottom$length==1)
  mtext(paste0(out.name, "_totij=", top.ij+bot.ij, "_topij=", top.ij, 
               "_botij=", bot.ij), 
        side=3, line=-1, adj=0.5, cex=0.5)
  dev.off()
  
  rm(df); gc()
  print(paste0(chr, " done!"), quote=FALSE)
}

# rm(list=ls())