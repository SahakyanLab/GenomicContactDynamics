################################################################################
# Plot fraction of long-range contacts per Cp per tissue. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/2_Expl_contacts"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# Count of long-range contacts
csv.dir = paste0(wk.dir, "/out_countLRcontacts")
out.dir = paste0(wk.dir, "/out_contactcount_plot")
### OTHER SETTINGS #############################################################
gcb.v = c("min2Mb", "min05Mb")
chr.v = paste0("chr", c(1:22, "X", "ALL")) 
Cp.v = 1:21
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(reshape2)
library(RColorBrewer)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
Cp.v <- as.character(sort(unique(Cp.v), decreasing=FALSE))
coul <- colorRampPalette(rev(brewer.pal(n=11,name="Spectral")))(length(Cp.v))
names(coul) <- Cp.v

for(gcb in gcb.v){
  
  print(paste0(gcb, "...."), quote=FALSE)
  
  for(chr in chr.v){
    
    LR <- read.csv(file=paste0(csv.dir, "/", gcb, "_", chr, "_HiCLRcontactsCount.csv"),
                   header=TRUE, row.names=1, stringsAsFactors=FALSE)
    LR <- data.matrix(LR)
    
    # Order ct from highest to lowest average Cp
    aveCp.v <- LR[rownames(LR)[rownames(LR)!="allCp"],]
    aveCp.v <- aveCp.v*as.numeric(rownames(aveCp.v))
    aveCp.v <- colSums(x=aveCp.v, na.rm=FALSE)
    
    tot.ij <- LR["allCp",names(aveCp.v)]
    if( length(aveCp.v)==length(tot.ij) ){
      ordr.ave <- sort(aveCp.v/tot.ij[names(aveCp.v)], decreasing=TRUE)
    } else {
      stop(paste0(gcb, " ", chr, ": Checkpoint 1."))
    }
    rm(aveCp.v)
  
    # Convert to fraction
    denom.mx <- matrix(data=LR["allCp",], nrow=nrow(LR)-1L, ncol=ncol(LR), byrow=TRUE)
    LR <- LR[rownames(LR)!="allCp",]/denom.mx
    rm(denom.mx)
    
    if( any(colSums(LR)!=1) ){
      stop(paste0(gcb, " ", chr, ": Fraction values don't add up to 1."))
    }
    
    # Order of ct from highest to lowest fr of Cp=21
    ordr.Cp21 <- sort(LR["21",], decreasing=TRUE)
    
    # Melt for plotting
    LR <- reshape2::melt(LR)
    colnames(LR) <- c("Cp", "ct", "fr")
    LR$Cp <- as.numeric(LR$Cp)
    LR$fr <- as.numeric(LR$fr)
    LR$ct <- as.character(LR$ct)
    LR$ct[LR$ct=="allCT"] <- "All"
    
    LR$Cp <- factor(LR$Cp, levels=Cp.v)
    
    # PLOT - Order ct from highest to lowest average Cp
    
    LR$ct <- factor( as.character(LR$ct), levels=unique(c(names(ordr.ave), "All")) )
     
    ggplot(data=LR, aes(x=ct, y=fr, fill=Cp)) +
      geom_col(position="fill") +
      scale_fill_manual(values=coul) +
      labs(title=paste0(gcb, "_", chr, "_LRcontacts_tissuehightolowWeightedAveCp"), 
           x=expression(bold("c"["p"])),
           y="Fraction of long-range contacts", fill="Cp") + 
      bgr2 +
      theme(axis.text.x=element_text(angle=90, size=18, hjust=1))

    ggsave(filename=paste0(out.dir, "/", gcb, "_", chr, "_LRcontacts_orderWAveCp_frplot.pdf"),
           unit="in", width=10, height=10)
    
    # PLOT - Order of ct from highest to lowest fr of Cp=21
    
    LR$ct <- factor( as.character(LR$ct), levels=unique(c(names(ordr.Cp21), "All")) ) 
    
    ggplot(data=LR, aes(x=ct, y=fr, fill=Cp)) +
      geom_col(position="fill") +
      scale_fill_manual(values=coul) +
      labs(title=paste0(gcb, "_", chr, "_LRcontacts_tissuehightolowfrCp21"), 
           x=expression(bold("c"["p"])),
           y="Fraction of long-range contacts", fill="Cp") + 
      bgr2 +
      theme(axis.text.x=element_text(angle=90, size=18, hjust=1))
    
    ggsave(filename=paste0(out.dir, "/", gcb, "_", chr, "_LRcontacts_orderfrCp21_frplot.pdf"),
           unit="in", width=10, height=10)
    
    rm(LR)
    gc()
    
    print(paste0(chr, " done!"), quote=FALSE)
    
  }
  
}

# rm(list=ls()); gc()
