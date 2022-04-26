################################################################################
# Visualise contact maps using different contact value metrics
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
    prefix = "/Users/ltamon"
    lib = paste0(prefix, "/DPhil/lib")
    data.dir = paste0(prefix, "/Database")
    wk.dir = paste0(prefix, "/SahakyanLab/GenomicContactDynamics/13_Simulation")
    CII.dir = paste0(prefix, "/SahakyanLab/GenomicContactDynamics/11_Complementarity/z_ignore_git")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    prefix = "/project/sahakyanlab/ltamon"
    lib = paste0(prefix, "/DPhil/lib")
    data.dir = paste0(prefix, "/Database")
    wk.dir = paste0(prefix, "/DPhil/GenomicContactDynamics/21_Simulation")
    CII.dir = paste0(prefix, "/DPhil/GenomicContactDynamics/11_Constraints")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
  
}

param.file = paste0(wk.dir, "/param.csv") #"/C2_generate_map/param.csv")
param.v <- read.csv(file=param.file, stringsAsFactors=F, header=T)
param.ind = 12
param.v <- param.v[param.ind,]

species.id = param.v[["species.id"]]

if(species.id=="human"){
  
  simmap.dir = paste0(wk.dir, "/sim_maps")
  Cs.raw.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
  Cs.norm.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
  CII.dir = paste0(CII.dir, "/out_constraints/merged_final")
  chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
  bin.len = 40000
  gcb = "min2Mb"
  chr.v = paste0("chr", c(1:22, "X"))
  
} else if(species.id=="ath"){
  
  Cs.norm.dir = paste0(data.dir, "/arabidopsis_HiC/Wang2015_20kb_combined_contacts/HiCNorm_Wang")
  CII.dir = paste0(CII.dir, "/out_constraints_ath_20kb/merged_final")
  chrlen.file = paste0(data.dir, "/genome_info/Ath_TAIR10_chr_info.txt") 
  
} else if(species.id=="osa"){
  
  Cs.norm.dir = paste0(data.dir, "/rice_HiC/Liu2017_50kb_combined_contacts/ICE_Liu")
  CII.dir = paste0(CII.dir, "/out_constraints_osa_50kb")
  chrlen.file = paste0(data.dir, "/genome_info/Osa_IRGSP1.0_chr_info.txt")
  
} else if(species.id=="dme"){
  
  Cs.norm.dir = paste0(data.dir, "/drosophila_Chathoth2019_HiC/Chathoth2019_10kb_combined_contacts/KR_Chathoth")
  CII.dir = paste0(CII.dir, "/out_constraints_dme_10kb")
  chrlen.file = paste0(data.dir, "/genome_info/Dme_dm6_chr_info.txt")
  bin.len = 10000
  gcb = "min0Mb"
  chr.v = paste0("chr", c("2L", "2R", "3L", "3R", "4", "X")) 
  
} else {
  stop("Species data not available.")
}

CII.disc.kmer.5.dir = CII.disc.align.5.dir = CII.disc.G.5.dir = CII.dir
CII.cont.kmer.5.dir = CII.cont.align.5.dir = CII.cont.G.5.dir = CII.dir
CII.disc.kmer.10.dir = CII.disc.align.15.dir = CII.disc.G.15.dir = CII.dir
CII.cont.kmer.15.dir = CII.cont.align.15.dir = CII.cont.G.15.dir = CII.dir
out.dir = paste0(wk.dir, "/out_generate_map_gcdmns")
### OTHER SETTINGS #############################################################
#gcb = "min2Mb" #"min0Mb" for ath
#bin.len = 40000 #40000 #20000 #50000 

#-------------------SELECT CONTACT MAPS

#chr.v = c("chr3R", "chr3L") #"chr3L" #paste0("chr", c("2L", "2R", "3L", "3R", "4", "X")) 

# Map id format: <cell/tissue>-<metric name>. Metric name should match source 
# directory name, e.g. for metric name Cs.norm directory is Cs.norm.dir. 
# For c||, <CII>.<disc/cont>.<kmer/align>.<(0,100)>, e.g. 5 in "CII.disc.kmer.5" 
# is the cutoff percentage for categorisation. disc means discrete (categorised CII), 
# cont means continuouos (orig CII).

# Specify metric for upper and lower matrix by writing element of ct.v and 
# metric.v as <cell type/metric upper>;<cell type/metric lower>. 

#ct.v = c("All;BG3", "All;BG3",
#         "All;Kc167", "All;Kc167") #param.v[["ct.v"]]
#metric.v = c("CII.cont.kmer.5;Cs.norm", "CII.disc.kmer.5;Cs.norm",
#             "CII.cont.kmer.5;Cs.norm", "CII.disc.kmer.5;Cs.norm") #param.v[["metric.v"]]

#ct.v = c("All;osa", "All;osa")
#metric.v = c("CII.cont.kmer.5;Cs.norm", "CII.disc.kmer.5;Cs.norm")

#ct.v = c("All;ath", "All;ath")
#metric.v = c("CII.cont.kmer.5;Cs.norm", "CII.disc.kmer.5;Cs.norm")

ct.v = param.v[["ct.v"]]
metric.v = param.v[["metric.v"]]
 
# Useful in case not whole chr is to be plotted
out.id = gsub(x=paste(paste(ct.v, metric.v, sep="_"), collapse="_"), 
              pattern=".", replacement="", fixed=T)
out.id = gsub(x=out.id, pattern=";", replacement="", fixed=T)
out.id = paste0(out.id, collapse="_")
chr.id = ifelse(length(chr.v)==1, chr.v, paste(chr.v[c(1, length(chr.v))], collapse="To"))
out.id = paste0(species.id, "_", chr.id, "_", out.id) #paste0(species.id, "_chr122X_", out.id) 

if( length(ct.v)!=length(metric.v) ){
  
  stop("Each element of ct.v and metric.v should correspond such that the two 
       vectors have the same length.")

} else {
  
  map.id.v <- expand.grid(chr.v, paste0(ct.v, "-", metric.v))
  map.id.v <- apply(X=map.id.v, MARGIN=1, FUN=paste, collapse="-")
  
}

scaleContactByDist.TF = param.v[["scaleContactByDist.TF"]]
scaled.disc.cutoff = as.numeric(param.v[["scaled.disc.cutoff"]])

# Convert metric values to contact probability?
contProb = F

#-------------------FILTER CONTACTS

# If both incl.bin.x and incl.bin.y lists are NULL, use whole chr.
# Upper triangle perspective, i -> y, j -> x
incl.bin.x = NULL
incl.bin.y = NULL
mask.bin.x = NULL #list(3563:6232) #NULL
mask.bin.y = NULL #list(1:3563) #NULL
# If closed vector gap.range is NULL, no filtering. 
gap.range = c(50, Inf)

#-------------------SET PLOT PARAMETERS

# Controls the x and y axes bounds; upper triangle perspective
limits.x = NULL 
limits.y = NULL 
# Plot values symmetrically e.g. plot value for (1,2) and (2,1)? Depends on limits set.
symmetric = T

# Mark bins along x- or/and y-axis
#tmp = seq(1000, 1500, 100)
mark.x = NULL #c(1, tmp, 2812-tmp+1, 2812)
mark.y = NULL #mark.x
#rm(tmp)

# Output specifications

# If scalebr.v==NULL, no scale bar
# scalebr.v = c(xmin=1, xmax=100, ymin=1, ymax=50)
scalebr.v = c(xmin=1, xmax=4e6/bin.len, ymin=1, ymax=50)
res = 100
# Number of rows and columns plot will be displayed
out.dim = c(nrow=1, ncol=2)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(compiler)
library(cowplot)
library(data.table)
library(ggnewscale)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(scales)
library(splines)
library(foreach)
library(itertools)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/simulation_lib/convertToContactProb.R"))
source(paste0(lib, "/simulation_lib/getmapdir.R"))
source(paste0(lib, "/simulation_lib/getContactDF.R"))
source(paste0(lib, "/simulation_lib/processForMap.R"))
source(paste0(lib, "/simulation_lib/makeMatrixMap.R"))
source(paste0(lib, "/simulation_lib/filterContacts.R"))
source(paste0(lib, "/simulation_lib/contactprobVsDistance.R"))
source(paste0(lib, "/categoriseValues.R"))
source(paste0(wk.dir, "/lib/scaleContactByDist.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( any(duplicated(map.id.v)) ){
  warning("Duplicated map in map.id.v.")
}

tmp <- list(gap.range, incl.bin.x, incl.bin.y, mask.bin.x, mask.bin.y)
            #limits.x, limits.y, mark.x, mark.y)
names(tmp) <- c("gap", "inclx", "incly", "maskx", "masky")
                #"limitsx", "limitsy", "markx", "marky")

for( nme in names(tmp) ){
  out.id <- paste0(out.id, "_", nme, paste(tmp[[nme]], collapse="To"))
}

out.name0 <- paste0(gcb, "_", out.id, "_res", res)
out.name0 <- paste0(out.name0, "_scaleContactByDist", scaleContactByDist.TF,
                    "_scaleddisc", scaled.disc.cutoff)

p.lst <- list()
len <- length(map.id.v)
for(M in 1:len){
  
  map.id <- map.id.v[M]
  chr <- strsplit(x=map.id, split="-", fixed=T)[[1]][1]
  ct.p <- strsplit(x=map.id, split="-", fixed=T)[[1]][2]
  metric.p <- strsplit(x=map.id, split="-", fixed=T)[[1]][3]
  
  map.v <- paste(strsplit(x=ct.p, split=";", fixed=T)[[1]],
                 strsplit(x=metric.p, split=";", fixed=T)[[1]], sep="-") 
  map.v.len <- length(unique(map.v))
  
  if( !length(map.v.len)%in%c(1,2) ){
    stop(paste0(map.id, ": Checkpoint 1."))
  }
  
  df <- list()
  for(m in 1:map.v.len){
    
    map <- map.v[m]
    ct <- strsplit(x=map, split="-", fixed=T)[[1]][1]
    metric <- strsplit(x=map, split="-", fixed=T)[[1]][2]
    metric.dir <- getmapdir(metric=metric, simmap.dir=simmap.dir)
    
    # Upper triangle contacts only, df$value can have 0
    df[[map]] <- getContactDF(metric.dir=metric.dir, metric=metric, 
                              gcb=gcb, chr=chr, ct=ct, gap.range=gap.range, 
                              incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y, 
                              mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y,
                              chrlen.file=chrlen.file, bin.len=bin.len, 
                              invalidij.action=NA, species.id=species.id)[,c("i", "j", "value")]
                            
    rownames(df[[map]]) <- NULL
    
    if(m==2){
      
      # Convert to lower triangle contacts
      df[[2]] <- df[[2]][,c("j", "i", "value")]
      colnames(df[[2]]) <- c("i", "j", "value")
      
    }
    
    # Convert to contact probability by dividing by max value per bin
    if( !grepl(x=metric, pattern="CII.disc.") & contProb){
      
      df[[map]] <- convertToContactProb(df=df[[map]], metric=metric, 
                                        tot.bin=length(unique(c(df[[map]]$i, df[[map]]$j)))
                                        )
      
    }
    
    # For metrics other than CII, 0 values converted to NAs
    # CII continuous can have 0 if sequences are identical but
    # but check whether this happens.
    
    if( !grepl(x=metric, pattern="CII.disc.|CII.cont.") ){
      df[[map]][ !is.na(df[[map]]$value) & df[[map]]$value==0,"value" ] <- NA
    } else if( grepl(x=metric, pattern="CII.cont.") ){
      
      if( max(df[[map]]$value, na.rm=T)>=0 ){
        stop(paste0(map, ": 0 and/or positive values in CII continuous."))
      }
      
    }
    
    if( !grepl(x=metric, pattern="CII.disc.") & any(na.omit(df[[map]]$value)==0) ){
      
      rm(df)
      stop(paste0(map, ": 0s present."))
      
    }
    
    rm(ct, metric, metric.dir)
    
  } # map.v.len for loop end
  
  if( symmetric & length(df)==1 ){
    
    # Use upper triangle to get lower triangle
    df[[2]] <- df[[1]][,c("j", "i", "value")]
    colnames(df[[2]]) <- c("i", "j", "value")

  }
  
  out.name <- paste(out.name0, paste0("contProb", contProb), chr, 
                    paste(map.v, collapse="_"), sep="_")
  
  if(scaleContactByDist.TF){
    
    df <- scaleContactByDist(df.lst=df, bin.len=bin.len, 
                             out.filepath=paste0(out.dir, "/", out.name, "_scaleContactByDistPlot"),
                             plot.title=paste0(out.name, "_invalidij.actionNA"))
    
    subj.ind <- which(grepl(x=names(df), pattern="CII.cont", fixed=T))
    if(scaled.disc.cutoff>0){
      
      df[[subj.ind]]$value <- categoriseValues(val.v=df[[subj.ind]]$value, cutoff=scaled.disc.cutoff)
      names(df)[subj.ind] <- "All-CII.disc.kmer.5"
      metric.p <- "CII.disc.kmer.5;Cs.norm"
      
    }
    
    if(chr=="chr1"){
      
      for( ind in 1:length(df) ){
        
        incl.TF <- filterContacts(ij.df=df[[ind]][,c("i","j")], gap.range=gap.range,
                                  incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y,  
                                  mask.bin.x=list(3563:6232), mask.bin.y=list(1:3563))
        
        df[[ind]]$value[!incl.TF] <- NA
        
      }
      
    }
      
  }

  # Plot
  p.lst[[paste0(M, map.id)]] <- makeMatrixMap(df.lst=df, check.dup=F, symmetric=symmetric,
                                              metric.v=strsplit(x=metric.p, split=";", fixed=T)[[1]][1:length(df)],
                                              plot.title=paste0(out.name, "_invalidij.actionNA"), 
                                              scalebr.v=scalebr.v, mark.x=mark.x, mark.y=mark.y,
                                              limits.x=limits.x, limits.y=limits.y, 
                                              is.contProb=contProb, is.scaleContactByDist=scaleContactByDist.TF,
                                              species.id=species.id)
                                   
  print(paste0(out.name, " done!"), quote=F)
  #rm(df, map.v, chr, ct.p, metric.p, out.name)
  gc()
  
} # map.id.v for loop end

p.lst.len <- length(p.lst)

foreach(inds=isplitVector(x=1:p.lst.len, chunkSize=prod(out.dim)), .inorder=T
        
) %do% {
  
  p.arr <- cowplot::plot_grid(plotlist=p.lst[inds], nrow=out.dim[1], ncol=out.dim[2],
                              align="none", axis="r", rel_widths=c(1,1), rel_heights=c(1,1),
                              labels=NULL, byrow=T)
  plot.ind <- ifelse(inds[1] < 10, paste0("0", inds[1]), inds[1])
  cowplot::save_plot(p.arr, filename=paste0(out.dir, "/", out.name0, "_contProb", contProb, 
                                            "_", plot.ind, ".png"), 
                     base_height=out.dim[1]*15, base_width=out.dim[2]*15, limitsize=F)
  
}

# rm(list=ls()); gc()

