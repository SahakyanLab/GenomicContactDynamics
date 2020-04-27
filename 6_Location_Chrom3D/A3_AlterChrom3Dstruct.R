################################################################################
# Customise the colour of the chromosomes in the model by altering the cmm file.
# Mac, R/3.6.1
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/6_Chrom3D"
    os = "Mac"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "IMR90_LMNB1_GSE49341_hg19" # "IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
cmm.dir = paste0(wk.dir, "/0_cmm")
out.dir = paste0(wk.dir, "/out_cmmModified")
### OTHER SETTINGS #############################################################
out.name = "rad0.03"
ploidy = "haploid"
iteration.v = c("3750", "48750", "90000", "375000", "750000", "1125000", 
                "1500000", "1875000", "2250000", "2625000", "3000000")
# If TRUE, modify radius of linker and marker (bead)
changeRadius = TRUE
# Two options to modify radius:
## If radiusConstant = TRUE, it will use the set linkerRadius and beadRadius below
## If radiusConstant = FALSE, it will use the set linkerRadius and beadRadius as
## multiplier to the original radii
radiusConstant = TRUE
linkerRadius = beadRadius = 0.03
# Per chromosome
changeColour <- TRUE
chr.order <- c(1:22,"X")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(RColorBrewer)
library(foreach)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
cmm.dir <- paste0(cmm.dir, "/", model.id, "/pm_", ploidy)
for(itr in iteration.v){
  
  # Load cmm file from Chrom3D containing the xyz coordinates and radius
  cmmfile <- paste0( cmm.dir, "/",
                     dir(cmm.dir, pattern=paste0(itr, ".cmm")) ) 
  cmm <- readLines(cmmfile)
  cmm.len <- length(cmm) 
  
  #--------------------RADIUS-------------------------------
  #--------------------RADIUS-------------------------------
  
  if(changeRadius){
    
    for( pat in c("<link id1", "<marker id") ){
      ind <- grep(pattern=pat, x=cmm)
      splt <- strsplit(x=cmm[ind], split='"')
      
      # Check if linker/bead radius similar to other values like that of colour
      if( length(unique( unlist( lapply(X=splt, FUN=length), use.names=FALSE ) ))!=1 ){
        stop("Checkpoint 1")
      }
      
      # Extract old radii 
      oldRadii <- as.numeric( 
        lapply(X=splt, FUN=function(lne){
          # Radius value
          lne[which(lne==" radius=") + 1]
        }) )
      
      rm("splt")
      
      ind.len <- length(ind)
      oldRadii.len <- length(oldRadii)
      
      ifelse(pat=="<link id1", newRadius <- linkerRadius, newRadius <- beadRadius)
      
      # Replace old radii
      if(ind.len==oldRadii.len){
        
        for(i in 1:ind.len){
          
          if(radiusConstant==FALSE){
            cmm[ ind[i] ] <- sub( pattern=oldRadii[i], 
                                  replacement=oldRadii[i]*as.numeric(newRadius),
                                  x=cmm[ ind[i] ] )
          } else {
            cmm[ ind[i] ] <- sub( pattern=oldRadii[i], 
                                  replacement=newRadius,
                                  x=cmm[ ind[i] ] )
          }
          
        }
        
      }
      
    }
    
    rm("oldRadii", "oldRadii.len", "ind", "ind.len"); gc()
    
  } # ChangeRadius end for loop

  #--------------------COLOURPERCHR-------------------------
  #--------------------COLOURPERCHR-------------------------
  
  if(changeColour==TRUE){
    
    # Indices of beads
    bead.ind <- grep(pattern="<marker id", x=cmm)
    bound.start <- c()
    counter <- 1L
    # Remove last cmm line </marker_set>
    bead.ind.len <- length(bead.ind)-1
    
    for(b in 1:bead.ind.len){
      if((bead.ind[b]+1)==bead.ind[b+1]){
        bound.start[counter] <- bead.ind[b]
        counter <- counter + 1L
      }
    }
    
    splt <- strsplit(x=cmm[bound.start], split='"')
    chr.v <- unlist( lapply(X=splt, FUN=function(lne){
      # chrID 
      lne[which(lne==" chrID=") + 1]
      chr.num <- gsub(pattern="chr|\\_|A|B", replacement="", x=lne[which(lne==" chrID=") + 1])
    }), use.names=FALSE
    )
    
    bound.end <- c(bound.start[-1]-1, cmm.len-1)
    mx <- cbind(start=bound.start, end=bound.end, chr=chr.v)
    mx <- mx[order(match(mx[,"chr"], chr.order)),]
    
    rm("bead.ind", "bead.ind.len", "bound.start", "bound.end")
    
    coul0 <- brewer.pal(11, "Spectral")
    coul <- colorRampPalette(rev(coul0))(length(chr.order))
    coul.float <- t( apply(X=col2rgb(coul), MARGIN=c(1,2), FUN=function(x){x/255}) )
    if(ploidy=="diploid"){
      coul.float <- coul.float[rep(1:nrow(coul.float), each=2),]
    }
    
    chr.v.len <- length(chr.v)
    splt.cmm <- strsplit(x=cmm, split='"')
    
    cmm.color.mod <- foreach(i=1:chr.v.len, .inorder=TRUE, .combine="c"
                             
    ) %do% {
      splt.cmm.sub <- splt.cmm[ mx[i,"start"]:mx[i,"end"] ]
      splt.cmm.sub.edit <- lapply(X=splt.cmm.sub, FUN=function(lne){
        lne[which(lne==" r=") + 1] <- coul.float[i, "red"]
        lne[which(lne==" g=") + 1] <- coul.float[i, "green"]
        lne[which(lne==" b=") + 1] <- coul.float[i, "blue"]
        paste(lne, collapse='"')
      })
      return(unlist(splt.cmm.sub.edit, use.names=FALSE))
    }
    
    cmm.new <- c(cmm[1], cmm.color.mod, cmm[cmm.len])
    
    rm("cmm.color.mod", "coul.float", "splt.cmm", "mx")
    
    if(cmm.len!=length(cmm.new)){
      stop("Lengths of original and modified cmms do not match.")
    }
    
  } # changeColour end for loop
  
  write(x=cmm.new, file=paste0(out.dir, "/", model.id, "_", out.name, "_", itr, 
                               "_", ploidy, ".cmm"))
  rm("cmm.new"); gc()
  
} # iteration sapply end

# rm(list=ls())