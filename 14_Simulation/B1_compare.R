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
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
Cs.raw.dir = paste0(data.dir, "/GSE87112/combined_contacts/RAW_primary_cohort")
Cs.norm.dir = Cp.dir = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_primary_cohort")
CII.disc.kmer.5.dir = CII.cont.kmer.5.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/polished/11_Constraints/out_group"
#CII.disc.kmer.5.dir = paste0(wk.dir, "/data/CII")
#SIM.disc.kmer.5.dir = paste0(wk.dir, "/sim_3.2")
SIM.4.2.kmer.5.dir = paste0(wk.dir, "/sim_4.2")
out.dir = paste0(wk.dir, "/out_compare")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr1"
ct = "FC"
bin.len = 40000

# Metric name should match source directory name.
# 5 in "CII.disc.kmer.5" is the cutoff percentage for categorisation. disc means
# discrete (categorised CII), cont means continuouos (orig CII). 
# <CII/SIM>.<disc/cont>.<kmer/align>.<(0,100)>
metric.v = c(subj="SIM.4.2.kmer.5", ref="Cs.norm")
#metric.v = c(subj="Cs.raw", ref="Cs.norm")

# compareContactMx() arguments
# If list incl.ind is NULL, use whole chr. 
incl.ind = list(1000:3200)
mask.ind = NULL
# If vector gap.range is NULL, no filtering. 
gap.range = c(50, 1900)
# Minimum value to be contact
c.offsubj.v = c(-0.1, 
                seq(from=0, to=0.005, by=0.0001),
                seq(from=0.01, to=0.06, by=0.01))
c.offref.v = c(-0.1, seq(from=0, to=0.5, by=0.01), 
               seq(from=0.6, to=2, by=0.1),
               seq(from=3, to=10, by=1))
nCPU = 7L

out.id = "only_gap50To1900_inclu1000To3200_r2" #"whole_gap50To1900_inclu1150To3120"

makeBoxplotValues = TRUE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(compiler)
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(wk.dir, "/lib/getContactDF.R"))
source(paste0(wk.dir, "/lib/compareContactMx.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
id <- paste(paste(names(metric.v), metric.v, sep="_"), collapse="_")
out.name <- paste(gcb, chr, ct, out.id, id, sep="_"); rm(id)
print(paste0(out.name, "..."), quote=FALSE)

genome <- read.table(file=chrLenfile, stringsAsFactors=FALSE, header=TRUE,
                     colClasses=c("character", "integer", "integer"))
mx.len <- ceiling(genome$length.bp[genome$chromosome==chr]/bin.len)
if(mx.len!=genome$bins.40kb[genome$chromosome==chr]){
  stop("mx.len wrong")
}; rm(genome)
temp <- expand.grid(1:mx.len, 1:mx.len)
temp <- temp[temp[,1]>temp[,2],]

MX <- list()
if(makeBoxplotValues){
  png(file=paste0(out.dir, "/", out.name, ".png"), width=3000, height=3000, res=300)
  par(mfrow=c(2,2))
}
for(m in c("subj", "ref")){
  
  metric <- metric.v[m]
  eval(parse(text=paste0(
    'metric.dir <- ', metric, '.dir'
  )))
   
  # Upper triangle
  incl.ind.v <- unique(unlist(incl.ind))
  df <- getContactDF(metric.dir=metric.dir, metric=metric.v[m], gcb=gcb, chr=chr, 
                     ct=ct, bins.i=incl.ind.v, bins.j=incl.ind.v, gap.range=gap.range)
  rm(incl.ind.v)
  if( nrow(df)!=((mx.len*mx.len)-mx.len)/2 ){
    stop("Wrong number of upper triangle contacts.") 
  }
  # Change value of unwanted contacts to NA
  df[is.infinite(df$value),"value"] <- NA
  # Reorder
  df <- df[order(df$i, df$j),]
  # Converted to lower matrix df
  colnames(df) <- c("j", "i", "value")
  
  MX[[m]] <- matrix(data=NA, nrow=mx.len, ncol=mx.len)
  if( !identical(as.numeric(df$i), as.numeric(temp$Var1)) | 
      !identical(as.numeric(df$j), as.numeric(temp$Var2)) ){
    stop("Order of df wrong.")
  }
 
  MX[[m]][ lower.tri(MX[[m]], diag=FALSE) ] <- as.numeric(df$value)
  # Fill upper triangle only
  MX[[m]] <- t(MX[[m]])
  
  rm(df); gc()
  
  if(makeBoxplotValues){
    boxplot(x=as.vector(MX[[m]]), outline=FALSE, ylab="Value",
            main=paste0(out.name, "_", metric.v[m]), cex.main=0.5)
    boxplot(x=as.vector(MX[[m]]), outline=TRUE, ylab="Value",
            main=paste0(out.name, "_", metric.v[m]), cex.main=0.5)
  }
  
}; rm(temp); gc()
if(makeBoxplotValues){ dev.off() }

COMPIJMX <- compareContactMx(MXsubj=MX$subj, MXref=MX$ref, 
                             c.offsubj.v, c.offref.v,  
                             incl.ind=incl.ind, mask.ind=mask.ind, 
                             nCPU=nCPU)
rm(MX); gc()
write.csv(COMPIJMX, file=paste0(out.dir, "/", out.name, ".csv"),
          row.names=FALSE)

end.time <- Sys.time()
end.time-start.time 

# rm(list=ls()); gc()