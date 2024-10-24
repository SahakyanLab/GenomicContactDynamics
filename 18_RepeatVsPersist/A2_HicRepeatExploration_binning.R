################################################################################
# HiCRepeatExploration after binning repeats in age ranking
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = paste0(home.dir, "/DPhil/lib")

rep.group = "fam" # "fam" | "subfam" 
minelm.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group)
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
transposon.dir = paste0(wk.dir, "/out_transposon")
out.dir = paste0(wk.dir, "/out_HicRepeatExploration")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
nCPU = 2L
bin.size = -1 # If -1, use number of elements in agerank
transposon = "yes" # yes - transposon only, no - not transposon only, all - all
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

#col.nme <- ifelse(rep.group=="fam", "repFamily", "repName")
agerank <- read.csv(file=paste0(agerank.dir, "/rep", rep.group, ".csv"),
                    header=TRUE, stringsAsFactors=FALSE)[,"repName"]#; rm(col.nme)

# Transposon-based filtering
if( transposon%in%c("yes", "no") ){
  
  print(paste0("Transposon-based filtering...", transposon, "..."), quote=FALSE)
  load(paste0(transposon.dir, "/hg19repeats_Transposon_yes_no.RData"))
  agerank <- intersect(agerank, TRANSPOSON[[transposon]][[rep.group]])
  print(length(agerank), quote=FALSE)
  rm(TRANSPOSON)
  
}

if(bin.size==-1){ bin.size <- length(agerank) }

out.dir <- paste0(out.dir, "/", rep.group, bin.size)
if( !dir.exists(out.dir) ){
  dir.create(out.dir)
}

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
colnames(temp) <- c(paste0("rep", rep.group), "repName")
write.csv(temp, file=paste0(agerank.dir, "/rep", rep.group, bin.size, ".csv"),
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
    
    MINELM.MX <- cbind(ntis=ntis, MINELM.MX)
    dimnames(MINELM.MX)[[2]] <- c("ntis", el.nme.v)
    rm(min.el, ij.len, ntis)
    save(MINELM.MX, file=paste0(out.dir, "/", chr, "_MinElm_", gcb, ".RData"))
    rm(MINELM.MX); gc()
    
    print(paste0(chr, " done!"), quote=FALSE)
    
  } #  itr for loop end
  
}
### END OF PARALLEL EXECUTION ###

# rm(list=ls()); gc()