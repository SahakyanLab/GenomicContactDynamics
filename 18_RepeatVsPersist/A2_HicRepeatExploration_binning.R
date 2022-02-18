################################################################################
# HiCRepeatExploration after binning repeats in age ranking
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/18_RepeatVsPersist"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "subfam" # "fam" | "subfam" 
minelm.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group)
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
out.dir = paste0(wk.dir, "/out_HicRepeatExploration")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
nCPU = 1L
bin.size = 24L
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chr.v.len <- length(chr.v)

out.dir <- paste0(out.dir, "/", rep.group, bin.size)
if( !dir.exists(out.dir) ){
  dir.create(out.dir)
}

col.nme <- ifelse(rep.group=="fam", "repFamily", "repName")
agerank <- read.csv(file=paste0(agerank.dir, "/rep", rep.group, ".csv"),
                    header=TRUE, stringsAsFactors=FALSE)[,col.nme]; rm(col.nme)

# Bin repeat elements
el.lst <- el.nme.v <- agerank
el.len <- length(el.nme.v)
if(bin.size<=0){ "bin.size argument should be >=1." }
if(bin.size>1){
  el.nme.v <- as.character( 1:ceiling(el.len/bin.size) ) 
}
x <- seq.int(from=1, to=el.len, by=bin.size)
y <- c(x[-1]-1, el.len)
el.lst <- mapply(x=x, y=y, FUN=function(x,y){return(el.lst[x:y])}, SIMPLIFY=FALSE)
names(el.lst) <- el.nme.v
rm(x, y, el.len)

# Save as binned repeats as agerank file 
temp <- lapply(X=el.lst, FUN=paste, collapse=";")
temp <- stack(temp)
colnames(temp) <- c("repSubfam", "repName")
write.csv(temp, file=paste0(agerank.dir, "/repsubfam", bin.size, ".csv"),
          row.names=FALSE)
rm(temp); gc()

toExport <- c("chr.v", "minelm.dir", "gcb", "el.nme.v", "el.lst", "out.dir")
#### PARALLEL EXECUTION #########
foreach(itr=isplitVector(1:chr.v.len, chunks=nCPU), 
        .inorder=FALSE, .export=toExport, .noexport=ls()[!ls()%in%toExport]
        
) %op% {
  
  for(i in itr){
    
    chr <- chr.v[i]
    load(paste0(minelm.dir, "/", chr, "_MinElm_", gcb, ".RData"))
    ntis <- MINELM.MX[,"ntis"]
    ij.len <- length(MINELM.MX[,1])
    min.el <- colnames(MINELM.MX)
    min.el <- min.el[min.el!="ntis"]
   
    MINELM.MX <- foreach(el.nme=el.nme.v, .inorder=TRUE, .combine="cbind"
                         
    ) %do% {
      
      el.v <- min.el[ min.el%in%el.lst[[el.nme]] ]
      if(length(el.v)==1){
        chunk <- MINELM.MX[,el.v] 
      } else if(length(el.v)>1){
        chunk <- rowSums(x=MINELM.MX[,el.v])
      } else if(length(el.v)==0){
        chunk <- rep(x=0, times=ij.len) 
      } else {
        stop("Checkpoint 1.")
      }
      return(chunk)
      
    } # el.nme.v foreach loop end
    
    dimnames(MINELM.MX)[[2]] <- el.nme.v
    MINELM.MX <- cbind(ntis=ntis, MINELM.MX)
    rm(min.el, ij.len, ntis)
    save(MINELM.MX, file=paste0(out.dir, "/", chr, "_MinElm_", gcb, ".RData"))
    rm(MINELM.MX); gc()
    
    print(paste0(chr, " done!"), quote=FALSE)
    
  } #  itr for loop end
  
}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()