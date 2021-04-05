################################################################################
# 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    binmx.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/binmx/out_bindata_1perc_HiCNorm"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
meanCp.id = "Cp1To21" # Which BINWMEANCP.DF to use
meanCp.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binWeightedMeanCp/", meanCp.id)
mutCp.dir = paste0(wk.dir, "/out_mutCalcPerCp")
out.dir = paste0(wk.dir, "/out_mutCalcVsMeanCp_bpsc")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
data.id = "donor_centric_PCAWG_sigEperc30_MMR12" # "CosmicNCV" | "donor_centric_PCAWG_as_is" | "donor_centric_PCAWG_sigEperc30_MMR12"
src.id = "Hg19" # "Hg19" | "hg38ToHg19"

mut.v = c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

calc.v = "Nmsitenorm" #c("Tmut", "Nmsite", "TmutDIVNmsite", "Nmsitenorm")
calc.id = "Nmsitenorm" # "Allcalc" | "Nmsitenorm"

wm = "wmeanCp" #"wmeanCp0" | "wmeanCp"

addjitter = FALSE # Add jitter points to boxplot?

SPLIT = "SLIDE"# BIN | SLIDE
# Maximum number of bins with the same mean Cp0 and mean Cp value are 617 and 117
# bins only, respectively. Both are less than 1% of bins. So setting to percPerBoxplot
# any value greater than 0.01 is acceptable.
percPerBp = NULL # percPerBp of data in each boxplot, only applicable to split="BIN"

# split = "SLIDE" parameters
d = 1 # Radius of sliding window
numSLwind=20 # Number of sliding (overlapping) windows
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/splitNumericVector.R"))
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/addwmeanCpToMUTCPDF.R"))
source(paste0(wk.dir, "/lib/aggregateDF.R"))
source(paste0(wk.dir, "/lib/identifyAltHyp.R"))
source(paste0(wk.dir, "/lib/doMannWhitneyAndScatter.R"))
source(paste0(wk.dir, "/lib/makebp.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
split.id <- ifelse(SPLIT=="BIN", paste0("BIN_percPerBp", percPerBp), 
                   paste0("SLIDE_d", d, "_numSLwind", numSLwind))
out.id <- paste0(gcb, "_", data.id, "_", src.id, "_meanCpFrom", meanCp.id, "_", calc.id, "_", wm, "_", split.id)

calc.v.len <- length(calc.v)
mut.v.len <- length(mut.v)
pdf(file=paste0(out.dir, "/", out.id, "_addjitter", addjitter, "_boxplots.pdf"), 
    height=calc.v.len*5, width=mut.v.len*5)
par(mfcol=c(calc.v.len, mut.v.len))

p.lst <- list()
for(mut in mut.v){
  
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=TRUE)

  # Load MUTCP.DF
  mutcp.id <- paste0(gcb, "_", data.id, "_", src.id, "_", mut.id)
  load(file=paste0(mutCp.dir, "/", mutcp.id, "_mutCalcPerCp.RData"))
  rm(mutcp.id)
  # Add mean Cp of bins to MUTCP.DF
  MUTCP.DF <- addwmeanCpToMUTCPDF(MUTCP.DF=MUTCP.DF, meanCp.dir=meanCp.dir, gcb=gcb)

  # Rows with finite weighted mean Cp
  eval(parse(text=paste0(
    wm, '.TF <- is.finite(MUTCP.DF[[wm]])'
  )))
  #wmeanCp0.TF <- is.finite(MUTCP.DF$wmeanCp0)
  #wmeanCp.TF <- is.finite(MUTCP.DF$wmeanCp)
  
  for(calc in calc.v){
    
    calc.TF <- is.finite(MUTCP.DF[[calc]])
      
    #for( wm in c("wmeanCp0", "wmeanCp") ){
      
      # Select only bins with finite weighted mean Cp and calculation value
      eval(parse(text=paste0(
        'bin.TF <-  calc.TF & ', wm, '.TF'
      )))
      x <- MUTCP.DF[bin.TF, c(calc, wm)]
      x <- x[order(x[[wm]], decreasing=FALSE),]
      
      # Total bins valid
      tot.bin <- length(x[[calc]])
      
      if(SPLIT=="BIN"){
        
        # Bin weigted mean Cp (non-overlapping intervals)
        f <- ceiling( seq_along(x[[calc]]) / (tot.bin*percPerBp/100) ) 
        temp <- stack(split( x=x[[calc]], f=f ))
        if( !identical(x[[calc]], temp$values) ){ stop("Checkpoint 2.") }
        x <- cbind.data.frame(x, ind=temp$ind)
        rm(temp, f); gc()
        
        # Min and max values of weighted mean Cp
        ext <- paste0("min=", min(x[[wm]], na.rm=TRUE), "_",
                        "max=", max(x[[wm]], na.rm=TRUE))
        
      } else if(SPLIT=="SLIDE"){
        
        # Split weigted mean Cp (slide, overlapping intervals)
        boundsSL <- c( min(x[[wm]]), ceiling(max(x[[wm]])) )
        slide <- splitNumericVector(x=x[[wm]], d=d, action="SLIDE", numSLwind=numSLwind,
                                    boundsSL=boundsSL)
        slide$midpoints <- NULL
        x <- lapply(X=slide, FUN=function(ind){
          x[[calc]][ind]
        })
        x <- stack(x); colnames(x)[colnames(x)=="values"] <- calc
        ext <- paste0("boundsSL_", paste(boundsSL, collapse="-"))
        
        rm(slide, boundsSL); gc()
        
      } else {
        stop("Invalid split argument. BIN or SLIDE only.")
      }
        
      # Number of bins per interval 
      binPerInt <- table(x$ind)[levels(x$ind)]
      binPerInt <- paste(paste(names(binPerInt), binPerInt, sep="="), collapse="_")
      binPerInt <- paste0("Totvalidbins=", tot.bin, "\n", binPerInt)
      
      # Make boxplot
      makebp(df=x, x="ind", y=calc, xlab=paste0(wm, "-bin"), ylab=calc, 
             addjitter=addjitter, plot.id=paste0(out.id, "_", mut.id, "_", binPerInt))
      
      #-------------------
      # Quantify p-value relative to bin=1 (lowest weighted mean Cp)
      
      x$ind <- as.numeric(x$ind)
      df <- aggregateDF(ind=x$ind, values=x[[calc]])
      id <- paste0(out.id, "_", mut.id, "_", calc, "_", ext)
      p.lst[[id]] <- doMannWhitneyAndScatter(x=x, df=df, calc=calc, plot.id=paste0(id, "_", binPerInt))$p
      
      rm(df, id, x, ext, binPerInt); gc()
      
    #} # wm for loop end
  
    rm(bin.TF, tot.bin)
    
  } # calc.v for loop ned
  
  rm(MUTCP.DF, mut.id); gc()
  print(paste0(mut, " done!"), quote=FALSE)
  
} # mut.v for loop end

dev.off()

if(calc.v.len>1){
  # For Allcalc
  p.arr <- ggarrange(plotlist=p.lst, nrow=mut.v.len, ncol=calc.v.len, legend=NULL)
  ggexport(p.arr, width=calc.v.len*10, height=mut.v.len*10,
           filename=paste0(out.dir, "/", out.id, "_scatplot.pdf"))
} else if(calc.v.len==1){
  # For single calc
  p.arr <- ggarrange(plotlist=p.lst, ncol=mut.v.len, nrow=calc.v.len, legend=NULL)
  ggexport(p.arr, height=calc.v.len*10, width=mut.v.len*10,
           filename=paste0(out.dir, "/", out.id, "_scatplot.pdf"))
} else {
  stop("Are you joking?!")
}

# rm(list=ls()); gc()





