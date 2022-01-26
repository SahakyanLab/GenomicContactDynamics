################################################################################
# Plot gene annotation-related lengths vs Cp 
# Consider longest transcript of gene per Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/16_GeneVsPersist"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
genelist.dir = out.dir = paste0(wk.dir, "/out_LTrUniqueIDPerCp")
data.dir = paste0(wk.dir, "/out_geneLength")
out.dir = paste0(wk.dir, "/out_geneLengthPlot")
### OTHER SETTINGS #############################################################
# Annotation file prefix
anno.nme = "hg19anno"
# Just NM (mRNA) because 
refseq = "ALL"
gcb.v = "min2Mb" #c("min2Mb", "min05Mb")
# Cp
nCPU = 2

plotOnly = FALSE
showOutliers = FALSE

lengths.v <- c("TRANSCRIPT.L", "INTRONS.dev.EXONS",
               "EXONS.L", "MEAN.EXON.L",
               "INTRONS.L", "MEAN.INTRON.L")
repFree = TRUE

repPercPlot = FALSE
#lengths.v <- c("REP.PERC.TR.L", "REP.PERC.EX.L", "REP.PERC.INT.L")

config <- list(
  # c(1-colour, 2-prefix, 3-abbreviation)
  TRANSCRIPT.L = c("#FDC776", "L", "tr"),
  INTRONS.dev.EXONS = c("#FDC776", NA, NA),
  EXONS.L = c("#e25d6c", "L", "exon"),
  MEAN.EXON.L = c("#e25d6c", "mean L", "ex"),
  INTRONS.L = c("#3288BD", "L", "int"),
  MEAN.INTRON.L = c("#3288BD", "mean L", "int"),
  REP.PERC.TR.L = c("#FDC776", "%R", "tr"),
  REP.PERC.EX.L = c("#e25d6c", "%R", "ex"),
  REP.PERC.INT.L = c("#3288BD", "%R", "int")
)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools) # isplitVector
library(reshape)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# ANNOLENGTH.DF
load(file=paste0(data.dir, "/", anno.nme, "_ALL_annoLengths.RData"))

if(repFree==FALSE){
  RF.ind <- grep(x=colnames(ANNOLENGTH.DF), pattern=".RF", fixed=TRUE)
  ANNOLENGTH.DF <- ANNOLENGTH.DF[,-(RF.ind)]
}

lengths.v.len <- length(lengths.v)
# Including the repeat-free lenghts
lengths.v1 <- colnames(ANNOLENGTH.DF)

for(gcb in gcb.v){
  
  filenme <- ifelse(repPercPlot==FALSE, 
                    paste0(out.dir, "/", anno.nme, "_", refseq, "_", gcb, "_repFree", repFree),
                    paste0(out.dir, "/", anno.nme, "_", refseq, "_", gcb, "_rep"))
  filenme <- ifelse(showOutliers==FALSE, paste0(filenme, "_combBP.pdf"), 
                    paste0(filenme, "_outl_combBP.pdf") )
  
  pdf(file=filenme, width=23, height=ifelse(repPercPlot==FALSE, 30, 20)) 
  
  # Outer margins c(bottom, left, top, right) # oma=c(5, 0, 2.5, 0)
  par(mfrow=c(ceiling(lengths.v.len/2), 2), oma=c(7, 0, 2.5, 0))
  
  if(plotOnly==FALSE){
    
    # File with unique IDs of longest transcript per cp
    LtrUniqueIDstringPerCp <- readLines(con=paste0(genelist.dir, "/", gcb, "_", refseq, 
                                                   "_LtrUniqueIDstringPerCpStr"))
    cp.header.ind <- grep(x=LtrUniqueIDstringPerCp, pattern="_cp", fixed=TRUE)
    # Isolate cp numbers (did this to be consistent with the order in the file)
    cp.headers <- as.numeric( gsub(x=LtrUniqueIDstringPerCp[cp.header.ind], 
                                   pattern="[^0-9]", replacement="") )
    cp.uniqueIDs <- LtrUniqueIDstringPerCp[cp.header.ind+1L]
    rm(LtrUniqueIDstringPerCp)
    cp.uniqueIDs.len <- length(cp.uniqueIDs)
  }
  
  id <- paste0(anno.nme, "_", refseq, "_", gcb)
  
  for(i in 1:lengths.v.len){ 
    
    len <- lengths.v[i]
    lens <- lengths.v1[ grep(pattern=len, x=lengths.v1, fixed=TRUE) ]
   
    if(plotOnly==FALSE){
      
      # len and len (repeat-free)
      lens.dta <- ANNOLENGTH.DF[,lens]
      
      toExport <- c("cp.uniqueIDs", "cp.headers", "lens.dta", "len")
      LENCP.DF <- foreach(itr=isplitVector(x=1:cp.uniqueIDs.len, chunks=nCPU),
                                   .inorder=TRUE, .combine="rbind",
                                   .export=toExport, 
                                   .noexport=ls()[!ls()%in%toExport]
      ) %op% {
        
        chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
          uniqueID <- as.numeric( strsplit(x=cp.uniqueIDs[i], split=";")[[1]] )
          if( !len%in%c("INTRONS.dev.EXONS", "REP.PERC.TR.L",
                        "REP.PERC.EX.L", "REP.PERC.INT.L") ){
            lens.dta <- as.matrix(lens.dta)[uniqueID,]/1000
          }
          rm(uniqueID); gc()
          return( cbind( cp=rep(cp.headers[i]), lens.dta ) )
        })
        do.call("rbind", chunk)
      }
      
      dimnames(LENCP.DF)[[2]] <- c("cp", "L", "L.RF")[1:ncol(LENCP.DF)]
      LENCP.DF <- as.data.frame(LENCP.DF)
      # For making the boxplot
      LENCP.DF$cp <- with(LENCP.DF, 
                                   reorder(LENCP.DF$cp,
                                           LENCP.DF$cp, na.rm=T))
      save(LENCP.DF, file=paste0(out.dir, "/", anno.nme, "_", 
                                          refseq, "_", gcb, "_", len, ".RData"))
      
      
      rm(lens.dta, lens); gc()
      
    } else {
      load(file=paste0(out.dir, "/", anno.nme, "_", refseq, "_", gcb, "_", 
                       len, ".RData"))
    }
    
    # MEAN.INTRON.L (INTRONS.L/N.INTRONS) has NaNs due to genes having no introns 
    # N.INTRONS=0
    LENCP.DF <- LENCP.DF[is.finite(LENCP.DF[,"L"]),]
    # Melt df for plotting
    
    # Calculate once; the same for all lengths
    if(i==1){ 
      N.cp <- tabulate(LENCP.DF$cp)
      cp.v <- as.numeric( levels(LENCP.DF$cp) )
      cp.v.len <- length(cp.v)
    }
    
    # ncol.k to differentiate between LENCP.DF with 2 or 3 columns depending
    # if there's a repeat-free version of the length
    ncol.k <- ncol(LENCP.DF)
    if(ncol.k==3){ LENCP.DF <- melt(data=LENCP.DF, id.vars="cp") }
    
    # c(1-colour, 2-prefix, 3-abbreviation)
    config.len.v <- config[[len]]
    
    #jpeg(filename=ifelse(showOutliers==FALSE, 
    #                     paste0(out.dir, "/", id, "_", len, "_combBP.jpeg"),
    #                     paste0(out.dir, "/", id, "_", len, "_outl_combBP.jpeg")),
    #     units="in", width=13, height=10, res=500)
  
    # Default = par(mar=c(5,4,4,2)+0.1); c(bottom, left, top, right)
    # Default = par(gp=c(3,1,0)); c(axis title, axis labels, axis line)
    par(mar=c(5, 10, 2, 2)+0.1, mgp=c(3, 2, 0))
    
    # Main plot
    myplot <- boxplot(formula=as.formula( ifelse(ncol.k==3, 
                                                 "value~variable*LENCP.DF$cp", 
                                                 "L~cp" ) ), 
                      data=LENCP.DF, boxwex=0.4, 
                      outline=ifelse(showOutliers==FALSE, FALSE, TRUE),
                      ylab="", main="", cex.axis=3, 
                      col=c(config.len.v[1], "gray91")[1:(ncol.k-1L)], xaxt="n")
    
    # Constructing x-axis 
    if(ncol.k==3){ # With repeat-free length
      # Line separating two boxplots per Cp
      for(i in seq(0.5 , cp.v.len*2+1, 2)){ abline(v=i, lty=1, col="gray80")}
      
      # X-axis labels
      # "1" "1" "2" "2" "3" "3"
      xlabs <- sapply( strsplit(x=myplot$names , split='\\.') , 
                       FUN=function(x) x[[2]] )
      # "1" "2" "3" "4"
      xlabs <- xlabs[seq(1, length(xlabs), 2)]
      
      # X coordinates for axis labels 
      xcoord <- seq(1.5, cp.v.len*2 , 2)
      
    } else { # Without repeat-free length
      xlabs <- xcoord <- cp.v
    }
    axis(side=1, at=xcoord, labels=xlabs, tick=TRUE, cex.axis=3)
    
    # X-axis title
    #mtext(side=1, text=expression("c"["p"]), line=3, cex=2)
    
    # Y-axis title
    prefix <- config.len.v[2]
    len.abbr <- config.len.v[3]
    if(len=="INTRONS.dev.EXONS"){
      
      mtext(side=2, text=bquote("L"^"int"/"L"^"ex"), line=4.7, cex=3)
        
    } else {
      
      if(grepl(pattern="REP", x=len, fixed=TRUE)){
        mtext(side=2, text=bquote(.(prefix)^.(len.abbr)), 
              line=4.7, cex=3)
      } else {
        mtext(side=2, text=bquote(.(prefix)^.(len.abbr)~", "%*%10^3~"nt"), 
              line=4.7, cex=3)
      }
      
    }
    
    # Add N per cp
    #mtext(at=xcoord, side=1, text=as.character(N.cp), cex=0.7, line=4)
    
    if(i==lengths.v.len & repPercPlot==FALSE){
      #legend("topright", legend=c("full", "rep-free"), 
      #       col=c(config.len.v[1] , "gray91"),
      #       pch=15, bty="o", bg="white", pt.cex=3, 
      #       cex=1.5, horiz=F, inset=c(0, -0.13), xpd=TRUE)
      
      legend("topleft", legend=c("", "full", "", "rep-free"), 
             col=c("#FDC776", "#e25d6c", "#3288BD", "gray91"), # config.len.v[1]
             pch=15, bty="o", bg="white", pt.cex=3, 
             cex=3, horiz=FALSE, inset=.05, xpd=TRUE)
    }

    #dev.off()
    
    print(len)
    
    rm(LENCP.DF); gc()
    
  } # lengths.v.len for loop end
  
  # Add title for the page
  mtext(side=3, text=ifelse(showOutliers==FALSE, paste0(id, "_outline=NA"), id),
        outer=TRUE, line=1, cex=0)
  mtext(side=1, text=expression("c"["p"]), outer=TRUE, line=5, cex=4)
  
  dev.off()
  
  if(plotOnly==FALSE){ rm(cp.uniqueIDs, cp.headers); gc() }
  print(paste0(gcb, " done!"))
  
} # gcb.v for loop end

# rm(list=ls()); gc()


