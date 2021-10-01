################################################################################
# Compare contact matrices
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
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
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
simmap.dir = paste0(wk.dir, "/sim_maps")
Cs.raw.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
Cs.norm.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
CII.disc.kmer.5.dir = CII.cont.kmer.5.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/polished/11_Constraints/out_group"
out.dir = paste0(wk.dir, "/out_compare")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr1"
ct = "Ao"
bin.len = 40000

# Metric name should match source directory name.
# For c||, <CII>.<disc/cont>.<kmer/align>.<(0,100)>
# 5 in "CII.disc.kmer.5" is the cutoff percentage for categorisation. disc means
# discrete (categorised CII), cont means continuouos (orig CII). 
metric.v = c(subj="SIM.int.set2.1.0", ref="Cs.norm")

# Filtering contacts
# If both incl.bin.x and incl.bin.y lists are NULL, use whole chr. 
incl.x = 'incl.bin.x = NULL'
incl.y = 'incl.bin.y = NULL'
mask.x = 'mask.bin.x = list(3039:6232)' #'mask.bin.x = list(3039:6232)'
mask.y = 'mask.bin.y = list(1:3580)'  #'mask.bin.y = list(1:3563)' 
# If vector gap.range is NULL, no filtering. 
gap.v = 'gap.range = c(50, Inf)'
out.id = "whole_maskMidSquare_gap50up"

# Minimum value to be contact, input as string so I can print it out
c.subj = 'c.offsubj.v = c( -0.00005, 0, seq(0.0009,0.003,0.00005), seq(0.0035,0.005,0.00025))'
c.ref  = 'c.offref.v = c( -0.05, seq(0,5,0.05) )'
#c.subj = 'c.offsubj.v = c( -0.0001, seq(0,0.004,0.0001), seq(0.005,0.01,0.001), seq(0.02,0.06, 0.01) )'
#c.ref  = 'c.offref.v = c( -0.02, seq(0,1,0.02), seq(1.2,1.8,0.2), seq(2,20,2) )'
nCPU = 3L #~8G per core if boxplotOnly = FALSE

boxplotOnly = FALSE
makeBoxplotValues = TRUE
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
source(paste0(wk.dir, "/lib/getContactDF.R"))
source(paste0(wk.dir, "/lib/compareContactMx.R"))
source(paste0(wk.dir, "/lib/filterContacts.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
id <- paste(paste(names(metric.v), metric.v, sep="_"), collapse="_")
out.name <- paste(gcb, chr, ct, out.id, id, sep="_"); rm(id)
print(paste0(out.name, "..."), quote=FALSE)

if(boxplotOnly){ print("Only making boxplot...", quote=FALSE) }

print(incl.x, quote=FALSE) 
print(incl.y, quote=FALSE) 
print(mask.x, quote=FALSE) 
print(mask.y, quote=FALSE) 
print(gap.v, quote=FALSE) 
print(c.subj, quote=FALSE)
print(c.ref, quote=FALSE)

eval(parse(text=incl.x))
eval(parse(text=incl.y))
eval(parse(text=mask.x))
eval(parse(text=mask.y))
eval(parse(text=gap.v))
eval(parse(text=c.subj))
eval(parse(text=c.ref))
rm(incl.x, incl.y, mask.x, mask.y, gap.v, c.subj, c.ref)

genome <- read.table(file=chrLenfile, stringsAsFactors=FALSE, header=TRUE,
                     colClasses=c("character", "integer", "integer"))
mx.len <- ceiling(genome$length.bp[genome$chromosome==chr]/bin.len)
if(mx.len!=genome$bins.40kb[genome$chromosome==chr]){ stop("mx.len wrong") } 
rm(genome)
temp <- expand.grid(1:mx.len, 1:mx.len)
temp <- temp[temp[,1]>temp[,2],]

MX <- list()
if(makeBoxplotValues){
  png(file=paste0(out.dir, "/", out.name, ".png"), width=6000, height=3000, res=300)
  par(mfrow=c(2,4))
}

for(m in c("subj", "ref")){
  
  metric <- unname(metric.v[m])
  
  # Define contact map directory, metric.dir
  if( grepl(x=metric, pattern="SIM.", fixed=TRUE) ){
    metric.dir <- paste0(simmap.dir, "/", metric)
  } else {
    eval(parse(text=paste0(
      'metric.dir <- ', metric, '.dir'
    )))
  }
  
  # Get contact values, upper triangle only
  df <- getContactDF(metric.dir=metric.dir, metric=metric, 
                     gcb=gcb, chr=chr, ct=ct, gap.range=gap.range, 
                     incl.bin.x=incl.bin.x, incl.bin.y=incl.bin.y, 
                     mask.bin.x=mask.bin.x, mask.bin.y=mask.bin.y)
  if( nrow(df)!=((mx.len*mx.len)-mx.len)/2 ){
    stop("Wrong number of upper triangle contacts.") 
  }
  if( any(is.na(df$value)) & !grepl(x=metric, pattern="CII.|SIM.") ){
    stop("NA values in df.")
  }
  
  # Change value of unwanted contacts to NA
  df$value[df$include==0] <- NA
  
  # Divide cs-based values by sd 
  if( metric%in%c("Cs.raw", "Cs.norm") ){
    
    df$value <- (df$value)/sd(x=df$value, na.rm=TRUE)
    print("Scaling Cs values by sd...", quote=FALSE)
    
  }
  
  # Reorder contacts to follow upper triangle order
  df <- df[order(df$i, df$j),]
  # Converted to lower matrix df
  colnames(df) <- c("j", "i", "value")
  df$value <- as.numeric(df$value)
  
  # Convert to matrix format
  MX[[m]] <- matrix(data=NA, nrow=mx.len, ncol=mx.len)
  if( !identical(as.numeric(df$i), as.numeric(temp$Var1)) | 
      !identical(as.numeric(df$j), as.numeric(temp$Var2)) ){
    stop("Order of df wrong.")
  }
  
  MX[[m]][ lower.tri(MX[[m]], diag=FALSE) ] <- df$value
  # Fill upper triangle only so transpose
  MX[[m]] <- t(MX[[m]])
  
  #-------------------Fix cut-off values and make boxplots
  
  # Store for boxplot
  v <- sort(df$value, decreasing=FALSE)
  rm(df); gc()
  
  # Use min and max to bound cut-off ranges.
  max.v <- tail(unique(v), n=2)
  min.v <- min(v[v!=0], na.rm=TRUE)
  
  if( (sum(v%in%max.v)!=2) & metric%in%c("Cs.norm", "Cs.raw") ){
    print("Cut-off max checkpoint.")
  }
  
  # Add cut-off values such that [0, max(value)]
  if(metric!="Cp"){
    
    eval(parse(text=paste0(
      "c.off", m, ".v <- c.off", m, ".v[c.off", m, ".v<max.v[1]]; 
      c.off", m, ".v <- sort(c(c.off", m, ".v, max.v)); 
      last3 <- tail(c.off", m, ".v, n=3);
      c.off", m, ".v <- sort(unique(c(c.off", m, ".v, seq(last3[1], last3[2], length.out=10))))"
    )))
    rm(last3)
    
  }
  
  # Make boxplot
  if(makeBoxplotValues){
    
    p.title <- paste0(out.name, "_", metric.v[m], "_min=", min.v, "_max=", max.v[2])
    
    if( !grepl(x=metric, pattern="CII.") ){
      v.nZero <- v[v!=0 & !is.na(v)]
      boxplot(x=v.nZero, outline=FALSE, ylab="Value",
              main=paste0(p.title, "_no0sNoOutliers"), cex.main=0.5)
      boxplot(x=v.nZero, outline=TRUE, ylab="Value",
              main=paste0(p.title, "_no0sWithOutliers"), cex.main=0.5)
    }
    boxplot(x=v, outline=FALSE, ylab="Value",
            main=paste0(p.title, "_AllNoOutliers"), cex.main=0.5)
    boxplot(x=v, outline=TRUE, ylab="Value",
            main=paste0(p.title, "_AllWithOutliers"), cex.main=0.5)
    
    rm(v.nZero, p.title)
    
  }
  
  rm(max.v, min.v, v)
  
}
rm(temp); gc()
if(makeBoxplotValues){ dev.off() }

# Compare matrices
if(boxplotOnly==FALSE){
  COMPIJMX <- compareContactMx(MXsubj=MX$subj, MXref=MX$ref, 
                               c.offsubj.v=c.offsubj.v, 
                               c.offref.v=c.offref.v,  
                               incl.bin.x=NULL, incl.bin.y=NULL,
                               mask.bin.x=NULL, mask.bin.y=NULL,
                               gap.range=NULL, nCPU=nCPU)
  write.csv(COMPIJMX, file=paste0(out.dir, "/", out.name, ".csv"),
            row.names=FALSE)
}
rm(MX); gc()

end.time <- Sys.time()
end.time-start.time 

# rm(list=ls()); gc()
