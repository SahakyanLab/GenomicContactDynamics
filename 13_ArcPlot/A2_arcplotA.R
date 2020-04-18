################################################################################
# Arcplot showing 2Mb and 0.5 Mb contacts.
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
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(wk.dir, "/out_arcplotA")
### OTHER SETTINGS #############################################################
chr.v = paste0("chr", c(22:1, "X"))
# topCp = 3 would mean the top 3 Cps so 19:21
topCP = 1
cp.v = 21:1
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(R4RNA)
library(RColorBrewer)
library(yarrr)
library(grDevices)
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
cp <- cp.v[1:topCP]
for(chr in chr.v){
  if(plotOnly==FALSE){
    df <- list()
    # Load PERSIST.MX 0.5Mb version
    load(file=paste0(persist.dir, "/", chr, "_Persist_min05Mb.RData"))
    incl.TF <- PERSIST.MX$ntis%in%cp
    if(sum(incl.TF)==0){
      print(paste0(chr, " skipped!"), quote=FALSE)
      next
    }
    df$min05Mb <- cbind(PERSIST.MX$hits[incl.TF, c("i", "j")],
                        length=as.integer(1))
    rm(incl.TF, PERSIST.MX); gc()
    df$min05Mb$value <- df$min05Mb$j-df$min05Mb$i
    # Order by decreasing gap distance, for the assignment of colour
    df$min05Mb <- df$min05Mb[order(df$min05Mb$value, decreasing=TRUE),]
    #-------------------
    # Add colours
    df$min2Mb <- df$min05Mb[df$min05Mb$value>50,]
    # Contacts shared with 2Mb version gray, contacts unique to 0.5Mb version red
    df$min05Mb <- colourByValue(as.helix(df$min05Mb), breaks=c(0, 50, max(df$min05Mb$value)),  
                                cols=c(adjustcolor("darkred", alpha.f=0.2), 
                                       "gray70"),
                                log=FALSE, include.lowest=TRUE, right=TRUE, get=TRUE)
    if( sum(is.na(df$min05Mb$col))!=0 ){
      stop("Checkpoint: Range do not cover all values.")
    }
    df$min2Mb <- colourByUnknottedGroups(as.helix(df$min2Mb), 
                                         colorRampPalette(yarrr::piratepal("basel"))(63),
                                         get=TRUE)
    if( sum(is.na(df$min2Mb$col))!=0 ){
      stop("Checkpoint: Colors not enough for unknotted groups.")
    }
    save(df, file=paste0(out.dir, "/", chr, "_topCP", topCP, "_LRmap.RData"))
  } else {
    load(file=paste0(out.dir, "/", chr, "_topCP", topCP, "_LRmap.RData"))
  }
  
  pdf(file=paste0(out.dir, "/", chr, "_topCP", topCP, "_LRmap.pdf"), width=30, height=30)
  plotDoubleHelix(as.helix(df$min2Mb), as.helix(df$min05Mb), shape="circle", 
                  line=TRUE, arrow=TRUE)
  mtext("> 2 Mb", side=3, line=-40, adj=0.05)
  ij2Mb <- nrow(df$min2Mb)
  mtext(paste0(chr, "_topCP=", topCP, "_2Mbij=", ij2Mb, "_restij=", nrow(df$min05Mb)-ij2Mb), 
        side=3, line=-40, adj=0.5)
  mtext("> 0.5 Mb", side=1, line=-40, adj=0.05)
  
  dev.off()
  
  rm(df); gc()
  print(paste0(chr, " done!"), quote=FALSE)
}

# rm(list=ls())