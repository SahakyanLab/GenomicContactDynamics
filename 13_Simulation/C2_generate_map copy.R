################################################################################
# Visualise contact maps using different contact value metrics
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    prefix = "/Users/ltamon"
    lib = paste0(prefix, "/DPhil/lib")
    data.dir = paste0(prefix, "/Database")
    wk.dir = paste0(prefix, "/DPhil/GenomicContactDynamics/21_Simulation")
    CII.dir = paste0(prefix, "/DPhil/GCD_polished/11_Complementarity/z_ignore_git/out_constraints/merged_final")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    prefix = "/t1-data/user/ltamon" # "/stopgap/sahakyanlab/ltamon"
    lib = paste0(prefix, "/DPhil/lib")
    data.dir = paste0(prefix, "/Database")
    wk.dir = paste0(prefix, "/DPhil/GenomicContactDynamics/21_Simulation")
    CII.dir = paste0(prefix, "/DPhil/GenomicContactDynamics/pending/11_Constraints/out_constraints/merged_final")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
  
}
simmap.dir = paste0(wk.dir, "/sim_maps")
Cs.raw.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
Cs.norm.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
CII.disc.kmer.5.dir = CII.disc.align.5.dir = CII.disc.G.5.dir = CII.dir
CII.cont.kmer.5.dir = CII.cont.align.5.dir = CII.cont.G.5.dir = CII.dir
out.dir = paste0(wk.dir, "/out_generate_map")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
bin.len = 40000

#-------------------SELECT CONTACT MAPS

chr.v = paste0("chr", c(1))

# Map id format: <cell/tissue>-<metric name>. Metric name should match source 
# directory name, e.g. for metric name Cs.norm directory is Cs.norm.dir. 
# For c||, <CII>.<disc/cont>.<kmer/align>.<(0,100)>, e.g. 5 in "CII.disc.kmer.5" 
# is the cutoff percentage for categorisation. disc means discrete (categorised CII), 
# cont means continuouos (orig CII).

# Specify metric for upper and lower matrix by writing element of ct.v and 
# metric.v as <cell type/metric upper>;<cell type/metric lower>. 
ct.v = c("hg19;FC", "hg19;PM",
         "hg19;hg19", "hg19;hg19", 
         "hg19;hg19", "hg19;hg19",
         "hg19;FC", "hg19;PM",
         "hg19;FC", "hg19;PM")
metric.v = c("Cp;Cs.norm", "Cp;Cs.norm",
             "CII.cont.kmer.5;Cp", "CII.cont.align.5;Cp", 
             "CII.disc.kmer.5;Cp", "CII.disc.align.5;Cp", 
             "CII.cont.kmer.5;Cs.norm", "CII.cont.kmer.5;Cs.norm",
             "CII.disc.kmer.5;Cs.norm", "CII.disc.kmer.5;Cs.norm")

ct.v = c("hg19;hg19", "hg19;hg19", "FC;hg19", "hg19;FC", "hg19;FC",
         "hg19;hg19", "hg19;hg19", "FC;hg19", "hg19;FC", "hg19;FC",
         "hg19;hg19", "hg19;hg19", "FC;hg19", "hg19;FC", "hg19;FC")
metric.v = c("CII.cont.kmer.5;Cp",  "CII.disc.kmer.5;Cp",  "Cs.norm;Cp", "CII.cont.kmer.5;Cs.norm",  "CII.disc.kmer.5;Cs.norm",
             "CII.cont.align.5;Cp", "CII.disc.align.5;Cp", "Cs.norm;Cp", "CII.cont.align.5;Cs.norm", "CII.disc.align.5;Cs.norm",
             "CII.cont.G.5;Cp",     "CII.disc.G.5;Cp",     "Cs.norm;Cp", "CII.cont.G.5;Cs.norm",     "CII.disc.G.5;Cs.norm")

#ct.v = c("hg19;hg19")
#metric.v = c("CII.disc.kmer.5;Cp")
  
if( length(ct.v)!=length(metric.v) ){
  
  stop("Each element of ct.v and metric.v should correspond such that the two 
       vectors have the same length.")

} else {
  
  map.id.v <- expand.grid(chr.v, paste0(ct.v, "-", metric.v))
  map.id.v <- apply(X=map.id.v, MARGIN=1, FUN=paste, collapse="-")
  
}

# Convert metric values to contact probability?
contProb = F

#-------------------FILTER CONTACTS

# If both incl.bin.x and incl.bin.y lists are NULL, use whole chr.
# Upper triangle perspective, i -> y, j -> x
incl.x = 'incl.bin.x = NULL'
incl.y = 'incl.bin.y = NULL'
mask.x = 'mask.bin.x = NULL' #'mask.bin.x = list(3038:6232)' #'mask.bin.x = list(3039:6232)'
mask.y = 'mask.bin.y = NULL' #'mask.bin.y = list(1:3565)'  #'mask.bin.y = list(1:3563)' 
# If closed vector gap.range is NULL, no filtering. 
gap.v = 'gap.range = c(50, Inf)'

#-------------------SET PLOT PARAMETERS

# Controls the x and y axes bounds; upper triangle perspective
limits.x = NULL #c(1965, 1974) #c(3090, 3099) #c(1915,1924) #c(3090, 3099) #c(1725,1734) 
limits.y = NULL #c(485, 494) #c(630, 639) #c(1860,1869) #c(630, 639) #c(1675,1684) 
# Plot values symmetrically e.g. plot value for (1,2) and (2,1)? Depends on limits set.
symmetric = T

# Mark bins along x- or/and y-axis
mark.x = NULL #c(485,494,1860,1869) #1965:1974 #1915:1924 #3090:3099 #1725:1734  
mark.y = NULL #c(1915,1924,1965,1974) #485:494 #1860:1869 #630:639 #1675:1684 

# Output specifications

# Useful in case not whole chr is to be plotted
out.id = "chr1_supp_gcd_figure" #"chr17_whole_50ToInf_FC_CpCs" #"chr17_50ToInf_x1915To1924_y1860To1869_gene1_FC_CpCs" #"chr1_whole_50ToInf_FC_CpCIICs_confimation" #"genes1" #"SGIP1" #"KRAS" #"whole_maskMidSquare_gap50up_maskx3038To6232y1To3565"
# If scalebr.v==NULL, no scale bar
scalebr.v = c(xmin=1, xmax=100, ymin=1, ymax=50)
res = 300
# Number of rows and columns plot will be displayed
out.dim = c(nrow=3, ncol=5)
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
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/simulation_lib/convertToContactProb.R"))
source(paste0(lib, "/simulation_lib/getmapdir.R"))
source(paste0(lib, "/simulation_lib/getContactDF.R"))
source(paste0(lib, "/simulation_lib/processForMap.R"))
source(paste0(lib, "/simulation_lib/makeMatrixMap.R"))
source(paste0(lib, "/simulation_lib/filterContacts.R"))
source(paste0(lib, "/categoriseValues.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( any(duplicated(map.id.v)) ){
  warning("Duplicated map in map.id.v.")
}

out.name0 <- paste(gcb, out.id, res, sep="_")
                   
plot.id <- paste(incl.x, incl.y, mask.x, mask.y, gap.v, sep=";")
plot.id <- paste0(plot.id, "_invalidij.actionNA")

print(incl.x, quote=F) 
print(incl.y, quote=F) 
print(mask.x, quote=F) 
print(mask.y, quote=F) 
print(gap.v, quote=F) 

eval(parse(text=incl.x))
eval(parse(text=incl.y))
eval(parse(text=mask.x))
eval(parse(text=mask.y))
eval(parse(text=gap.v))

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
                              invalidij.action=NA)[,c("i", "j", "value")]
                            
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
    
    rm(ct, metric, metric.dir)
    
  } # map.v.len for loop end
  
  if( symmetric & length(df)==1 ){
    
    # Use upper triangle to get lower triangle
    df[[2]] <- df[[1]][,c("j", "i", "value")]
    colnames(df[[2]]) <- c("i", "j", "value")
    
  }
  out.name <- paste(out.name0, paste0("contProb", contProb), chr, 
                    paste(map.v, collapse="_"), sep="_")
  
  # Plot
  p.lst[[paste0(M, map.id)]] <- makeMatrixMap(df.lst=df, check.dup=F, symmetric=symmetric,
                                              metric.v=strsplit(x=metric.p, split=";", fixed=T)[[1]][1:length(df)],
                                              plot.title=paste0(out.name, "_", plot.id,
                                                                "_scale=", unname(scalebr.v["xmax"]-scalebr.v["xmin"])
                                              ), 
                                              scalebr.v=scalebr.v, mark.x=mark.x, mark.y=mark.y,
                                              limits.x=limits.x, limits.y=limits.y, contProb=contProb)
                                   
  print(paste0(out.name, " done!"), quote=F)
  rm(df, map.v, chr, ct.p, metric.p, out.name)
  gc()
  
} # map.id.v for loop end

p.arr <- cowplot::plot_grid(plotlist=p.lst, nrow=out.dim[1], ncol=out.dim[2],
                            align="none", axis="r", rel_widths=c(1,1), rel_heights=c(1,1),
                            labels=NULL, byrow=T)
cowplot::save_plot(p.arr, filename=paste0(out.dir, "/", out.name0, "_contProb", contProb, ".png"), 
                   base_height=out.dim[1]*15, base_width=out.dim[2]*15, limitsize=F)

# rm(list=ls()); gc()

