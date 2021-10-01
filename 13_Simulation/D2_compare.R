################################################################################
# Compare contact matrices
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
simmap.dir = paste0(wk.dir, "/sim_maps")
Cs.raw.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
Cs.norm.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
CII.disc.kmer.5.dir = CII.cont.kmer.5.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/polished/11_Constraints/out_group"
out.dir = paste0(wk.dir, "/out_boxplot")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr21"
bin.len = 40000

# Select contact maps

# Map id format: <cell/tissue>-<metric name>. Metric name should match source 
# directory name, e.g. for metric name Cs.norm directory is Cs.norm.dir. 
# For c||, <CII>.<disc/cont>.<kmer/align>.<(0,100)>, e.g. 5 in "CII.disc.kmer.5" 
# is the cutoff percentage for categorisation. disc means discrete (categorised CII), 
# cont means continuouos (orig CII). 
map.id.v = c(subj="Ao-Cs.norm", ref="AG-Cs.norm") # c(subj="SUBJREPLACE", ref="REFREPLACE")

# Filter contacts

# If both incl.bin.x and incl.bin.y lists are NULL, use whole chr. 
# Upper triangle perspective, i -> y, j -> x
incl.x = 'incl.bin.x = NULL'
incl.y = 'incl.bin.y = NULL'
mask.x = 'mask.bin.x = list(3038:6232)' #'mask.bin.x = list(3039:6232)'
mask.y = 'mask.bin.y = list(1:3565)'  #'mask.bin.y = list(1:3563)' 
# If vector gap.range is NULL, no filtering. 
gap.v = 'gap.range = c(50, Inf)'

out.id = "whole_maskMidSquare_gap50up_maskx3038To6232y1To3565_xlim0T06"

# Minimum value to be contact, input as string so I can print it out
c.subj = 'c.offsubj.v = c( -0.0001, seq(0,0.004,0.0001), seq(0.005,0.01,0.001), seq(0.02,0.1, 0.01) )'
c.ref  = 'c.offref.v = c( -0.05, seq(0,5,0.05) )'
#c.subj = 'c.offsubj.v = c( -0.0001, seq(0,0.004,0.0001), seq(0.005,0.01,0.001), seq(0.02,0.06, 0.01) )'
#c.ref  = 'c.offref.v = c( -0.02, seq(0,1,0.02), seq(1.2,1.8,0.2), seq(2,20,2) )'
nCPU = 3L # No. of cut=off combination
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(reshape2)
library(compiler)
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(wk.dir, "/lib/getmapdir.R"))
source(paste0(wk.dir, "/lib/getContactDF.R"))
source(paste0(wk.dir, "/lib/compareContactMx.R"))
source(paste0(wk.dir, "/lib/filterContacts.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
id <- paste(paste(names(map.id.v), map.id.v, sep="_"), collapse="_")
out.name <- paste(gcb, chr, out.id, id, sep="_")
rm(id)
print(paste0(out.name, "..."), quote=F)

print(incl.x, quote=F) 
print(incl.y, quote=F) 
print(mask.x, quote=F) 
print(mask.y, quote=F) 
print(gap.v, quote=F) 
print(c.subj, quote=F)
print(c.ref, quote=F)

eval(parse(text=incl.x))
eval(parse(text=incl.y))
eval(parse(text=mask.x))
eval(parse(text=mask.y))
eval(parse(text=gap.v))
eval(parse(text=c.subj))
eval(parse(text=c.ref))
rm(incl.x, incl.y, mask.x, mask.y, gap.v, c.subj, c.ref)

# Get total bins of chr
genome <- read.table(file=chrlen.file, stringsAsFactors=F, header=T,
                     colClasses=c("character", "integer", "integer"))
tot.bin <- ceiling(genome$length.bp[genome$chromosome==chr]/bin.len)
if(tot.bin!=genome$bins.40kb[genome$chromosome==chr]){
  stop("tot.bin wrong")
} 
rm(genome)

# Template for checking order of df
temp <- expand.grid(i=1:tot.bin, j=1:tot.bin)
temp <- temp[temp$i<temp$j,]
temp <- temp[order(temp$j, temp$i),]

MX <- list()
for(m in c("subj", "ref")){
  
  map.id <- map.id.v[m]
  ct <- strsplit(x=map.id, split="-", fixed=T)[[1]][1]
  metric <- strsplit(x=map.id, split="-", fixed=T)[[1]][2]
  
  metric.dir <- getmapdir(metric=metric, simmap.dir=simmap.dir)

  # Upper triangle contacts only
  df <- getContactDF(metric.dir=metric.dir, metric=metric, 
                     gcb=gcb, chr=chr, ct=ct, gap.range=gap.range, 
                     incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y, 
                     mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y,
                     chrlen.file=chrlen.file, bin.len=bin.len, invalidij.action=NA)

  # Convert df to matrix format
  df <- df[order(df$j, df$i, decreasing=FALSE),]
  
  if( !identical(as.numeric(df$i), as.numeric(temp$i)) | 
      !identical(as.numeric(df$j), as.numeric(temp$j)) ){
    stop("Order of df wrong.")
  }
  
  # Convert to matrix format; see ./lib/UpdfToUpTriMx.R for details 
  MX <- matrix(data=NA, nrow=tot.bin, ncol=tot.bin)
  MX[ upper.tri(MX, diag=F) ] <- df$value
  
  #-------------------Fix cut-off values 
  
  v <- sort(df$value, decreasing=F)
  rm(df); gc()
  
  # Use min and max to bound cut-off ranges.
  max.v <- tail(unique(v), n=2)
  if( (sum(v%in%max.v)!=2) & metric%in%c("Cs.norm", "Cs.raw") ){
    print("Cut-off max checkpoint.")
  }
  
  # Add cut-off values such that [0, max(value)]; 0 is not minimum if metric=CII
  if(metric!="Cp"){
    
    eval(parse(text=paste0(
      "c.off", m, ".v <- c.off", m, ".v[c.off", m, ".v<max.v[1]]; 
      c.off", m, ".v <- sort(c(c.off", m, ".v, max.v)); 
      last3 <- tail(c.off", m, ".v, n=3);
      c.off", m, ".v <- sort(unique(c(c.off", m, ".v, seq(last3[1], last3[2], length.out=10))))"
    )))
    rm(last3)
    
  }
  
  rm(max.v, v, map.id, metric, ct, metric.dir)
  
}

rm(temp); gc()

# Compare matrices
COMPIJMX <- compareContactMx(MXsubj=MX$subj, MXref=MX$ref, 
                             c.offsubj.v=c.offsubj.v, 
                             c.offref.v=c.offref.v,  
                             incl.bin.x=NULL, incl.bin.y=NULL,
                             mask.bin.x=NULL, mask.bin.y=NULL,
                             gap.range=NULL, nCPU=nCPU)
write.csv(COMPIJMX, file=paste0(out.dir, "/", out.name, ".csv"), row.names=F)
rm(MX); gc()

end.time <- Sys.time()
end.time-start.time 

# rm(list=ls()); gc()
