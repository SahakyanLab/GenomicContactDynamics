################################################################################
# Visualise contact maps using different contact value metrics
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
    CII.dir = "/Users/ltamon/DPhil/GCD_polished/11_Complementarity/out_group"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    #lib = "/t1-data/user/ltamon/DPhil/lib"
    #data.dir = "/t1-data/user/ltamon/Database"
    #wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
    lib = "/stopgap/sahakyanlab/ltamon/DPhil/lib"
    data.dir = "/stopgap/sahakyanlab/ltamon/Database"
    wk.dir = "/stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
    CII.dir = "/stopgap/sahakyanlab/ltamon/DPhil/GenomicContactDynamics/pending/11_Constraints/out_group"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    data.dir = "/home/ltamon/Database"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
  
}
simmap.dir = paste0(wk.dir, "/sim_maps")
Cs.raw.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
Cs.norm.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
CII.disc.kmer.5.dir = CII.cont.kmer.5.dir = CII.dir
out.dir = paste0(wk.dir, "/out_generate_map")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################

## Select metrics

#ct.v = "hg19" #"CTREPLACE"
#metric.v = sort(list.files(path=simmap.dir), decreasing=F)[-(1:3)]
#ct.v = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
#              "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC"))[1:3]
ct.v = "Bl" #c("FC", "FC", "FC")
metric.v = "CII.disc.kmer.5" #"CII.cont.kmer.5" #c("Cs.raw", "Cs.norm", "Cp")
# Map id format: <cell/tissue>-<metric name>. Metric name should match source 
# directory name, e.g. for metric name Cs.norm directory is Cs.norm.dir. 
# For c||, <CII>.<disc/cont>.<kmer/align>.<(0,100)>, e.g. 5 in "CII.disc.kmer.5" 
# is the cutoff percentage for categorisation. disc means discrete (categorised CII), 
# cont means continuouos (orig CII). 
map.id.v = paste(ct.v, metric.v, sep="-")

# Convert metric values to contact probability?
contProb = FALSE

## Filter contacts

gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
bin.len = 40000

# If both incl.bin.x and incl.bin.y lists are NULL, use whole chr.
# Upper triangle perspective, i -> y, j -> x
incl.x = 'incl.bin.x = NULL'
incl.y = 'incl.bin.y = NULL'
mask.x = 'mask.bin.x = NULL' #'mask.bin.x = list(3038:6232)' #'mask.bin.x = list(3039:6232)'
mask.y = 'mask.bin.y = NULL' #'mask.bin.y = list(1:3565)'  #'mask.bin.y = list(1:3563)' 
# If vector gap.range is NULL, no filtering. 
gap.v = 'gap.range = c(50, Inf)'
## Define what to do with metric value of invalid contacts. Useful for controlling
## what the values the legends of the heatmap will show.
## NA - set to NA; "drop" - remove those contacts; "none" - don't touch
#invalidij.action = NA

# Controls whether you'll get a square or symmetric map; upper triangle perspective
limits.x = NULL #c(1965, 1974) #c(3090, 3099) #c(1915,1924) #c(3090, 3099) #c(1725,1734) 
limits.y = NULL #c(485, 494) #c(630, 639) #c(1860,1869) #c(630, 639) #c(1675,1684) 
format = "symmetric" # symmetric | square

# Mark bins along x- or/and y-axis
mark.x = NULL #1965:1974 #1915:1924 #3090:3099 #1725:1734  
mark.y = NULL #485:494 #1860:1869 #630:639 #1675:1684 

## Output specifications

# Useful in case not whole chr is to be plotted
out.id = "wholeRdYlBu" #"genes1" #"SGIP1" #"KRAS" #"whole_maskMidSquare_gap50up_maskx3038To6232y1To3565"
# If scalebr.v==NULL, no scale bar
scalebr.v = NULL #c(xmin=830, xmax=839, ymin=996, ymax=1000) #c(xmin=1, xmax=100, ymin=6183, ymax=6232)
res = 300
# Number of rows and columns plot will be displayed
out.dim = c(nrow=6, ncol=4)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(compiler)
library(data.table)
library(reshape2)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/getmapdir.R"))
source(paste0(wk.dir, "/lib/getContactDF.R"))
source(paste0(wk.dir, "/lib/processForMap.R"))
source(paste0(wk.dir, "/lib/makeMatrixMap.R"))
source(paste0(wk.dir, "/lib/filterContacts.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name0 <- paste(gcb, out.id, format, res, paste(map.id.v, collapse="_"), sep="_")
                   
plot.id <- paste(incl.x, incl.y, mask.x, mask.y, gap.v, sep=";")
plot.id <- paste0(plot.id, "_invalidij.actionNA")

print(incl.x, quote=FALSE) 
print(incl.y, quote=FALSE) 
print(mask.x, quote=FALSE) 
print(mask.y, quote=FALSE) 
print(gap.v, quote=FALSE) 

eval(parse(text=incl.x))
eval(parse(text=incl.y))
eval(parse(text=mask.x))
eval(parse(text=mask.y))
eval(parse(text=gap.v))

#chrlen.df <- read.table(file=chrlen.file, stringsAsFactors=F, header=T,
#                        colClasses=c("character", "integer", "integer"))

p.lst <- list()
len <- length(map.id.v)
for(chr in chr.v){
  
  ## Get total bins of chr
  #tot.bin <- ceiling(chrlen.df$length.bp[chrlen.df$chromosome==chr]/bin.len)
  #if(tot.bin!=chrlen.df$bins.40kb[chrlen.df$chromosome==chr]){
  #  stop("tot.bin wrong")
  #} 

  ## Template for checking order of df
  #temp <- expand.grid(i=1:tot.bin, j=1:tot.bin)
  #temp <- temp[temp$i<temp$j,]
  #temp <- temp[order(temp$j, temp$i),]
  
  for(map.id in map.id.v){
    
    ct <- strsplit(x=map.id, split="-", fixed=T)[[1]][1]
    metric <- strsplit(x=map.id, split="-", fixed=T)[[1]][2]
    metric.dir <- getmapdir(metric=metric, simmap.dir=simmap.dir)
    
    out.name <- paste(out.name0, 
                      paste0("contProb", 
                             # Never use contact probability for CII discrete
                             ifelse(grepl(x=metric, pattern="CII.disc.", fixed=TRUE), FALSE, contProb)),
                      chr, ct, metric, sep="_")
    print(paste0(out.name, "..."), quote=F)
    
    # Upper triangle contacts only, df$value can have 0
    df <- getContactDF(metric.dir=metric.dir, metric=metric, 
                       gcb=gcb, chr=chr, ct=ct, gap.range=gap.range, 
                       incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y, 
                       mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y,
                       chrlen.file=chrlen.file, bin.len=bin.len, 
                       invalidij.action=NA)
    
    df$include <- NULL
    
    #-------------------Transform metric values
    
    # For metrics other than CII.disc., 0 values converted to NAs
    # CII continuous can have 0 if sequences are identical but
    # but check whether this happens.
    
    nonNA.TF <- !is.na(df$value)
    
    if( !grepl(x=metric, pattern="CII.disc.|CII.cont.") ){
      df[nonNA.TF & df$value==0,"value"] <- NA
    } else if( grepl(x=metric, pattern="CII.cont.") ){
      
      if( max(df$value, na.rm=TRUE)==0 ){
        stop(paste0(out.name, ": 0 in CII continuous."))
      }
      
    }
  
    # Convert to contact probability by dividing by max value per bin
    if( !grepl(x=metric, pattern="CII.disc.") & 
        #identical(as.numeric(df$i), as.numeric(temp$i)) & 
        #identical(as.numeric(df$j), as.numeric(temp$j)) & 
        contProb ){
      
      # Transform CII continuous to positive values that can be transformed
      # to probabilities. Did eumiro's answer so distribution will only be translated:
      # https://stackoverflow.com/questions/3931419/turn-a-negative-number-into-a-positive-for-probability
      if( grepl(x=metric, pattern="CII.cont.") ){
        
        # +1 so probability for minimum value will be 1 (not 0)
        df$value[nonNA.TF] <- df$value[nonNA.TF]-min(df$value, na.rm=TRUE)+1
        
      }
      
      # Convert df to matrix
      tot.bin <- length(unique(c(df$i, df$j)))
      MX <- matrix(data=NA, nrow=tot.bin, ncol=tot.bin)
      
      df <- df[order(df$i, df$j, decreasing=FALSE),]
      MX[ lower.tri(MX, diag=FALSE) ] <- df$value
      
      df <- df[order(df$j, df$i, decreasing=FALSE),]
      MX[ upper.tri(MX, diag=FALSE) ] <- df$value
      
      if( !isSymmetric(MX) ){
        stop(paste0(out.name, ": Matrix not symmetrical."))
      }
      
      # Rows that are not all NAs
      nonNArw.TF <- !apply( X=MX, MARGIN=1, FUN=function(rw) all(is.na(rw)) )
      
      # Maximum value per row (per bin)
      maxrw.v <- rep(NA, times=length(nonNArw.TF))
      maxrw.v[nonNArw.TF] <- apply(X=MX[nonNArw.TF,], MARGIN=1, FUN=max, na.rm=TRUE) 
      
      if( !identical(nonNArw.TF, !is.na(maxrw.v)) ){
        stop(paste0(out.name, ": Checkpoint 1."))
      }
      
      MX <- MX/maxrw.v
      df$value <- MX[ upper.tri(MX, diag=FALSE) ]
      
      rm(nonNArw.TF, MX, tot.bin, maxrw.v)
      
    }
    
    #-------------------Plot
    
    p.lst[[out.name]] <- makeMatrixMap(df=df, format=format, check.dup=FALSE, metric=metric,
                                       plot.title=paste0(out.name, "\n", plot.id,
                                                         "\nscale=", unname(scalebr.v["xmax"]-scalebr.v["xmin"])
                                       ), 
                                       scalebr.v=scalebr.v, mark.x=mark.x, mark.y=mark.y,
                                       limits.x=limits.x, limits.y=limits.y, contProb=contProb)
    
    print(paste0(out.name, " done!"), quote=FALSE)
    rm(df, metric.dir, ct, metric, out.name, nonNA.TF)
    gc()
    
  } # map.id.v for loop end
  
  #rm(tot.bin, temp)
  
} # chr.v for loop end

p.arr <- ggarrange(plotlist=p.lst, nrow=out.dim["nrow"], ncol=out.dim["ncol"], legend=NULL)
ggexport(p.arr, height=5*out.dim["nrow"]*res, width=10*out.dim["ncol"]*res, res=res,
         filename=paste0(out.dir, "/", out.name0, "_contProb", contProb, ".png"))
#ggexport(p.arr, height=10, width=20*m.v.len, filename=paste0(out.dir, "/", out.name, ".pdf"))

# rm(list=ls()); gc()

## Transformation to reduce positive skewness even for zero (was
## already converted to NA before transformation) and negative values
#if( grepl(x=metric, pattern="Cs.raw|Cs.norm|CII.cont.") ){

#  nonNA.TF <- !is.na(df$value)
#  df$value[nonNA.TF] <- log(df$value[nonNA.TF] + 1)
#  #df$value[nonNA.TF] <- (df$value[nonNA.TF])^(1/3)
#  rm(nonNA.TF)

#}
