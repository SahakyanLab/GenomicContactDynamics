################################################################################
# Identify potential unmappable areas (and the largest of which), areas that 
# don't form contacts in all tissues.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(wk.dir, "/out_reduceContactMxGap")
out.dir = paste0(wk.dir, "/out_unmappable")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
metric = "Cs.norm"
out.id = "whole"
ct.v = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
              "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC"))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
ct.v.len <- length(ct.v)
out.name <- paste(gcb, out.id, metric, sep="_")

unmap <- list()
unmapbig <- list()
for(chr in chr.v){
  
  data.name <- paste(gcb, chr, out.id, metric, sep="_")
  x <- read.csv(file=paste0(data.dir, "/", data.name, "_contactMxGap_reduced.csv"), 
                row.names="X")
  if( any(duplicated(x)) ){
    stop(paste0(chr, ": Duplicated rows."))
  }
  
  x <- x[order(x$red.width, x$ct, decreasing=T),]
 
  # Unmappble area id
  sort.id <- paste0(x$red.start, ":", x$red.end)
  
  # Width per unique sort.id
  width.v <- x$red.width
  names(width.v) <- sort.id
  width.v <- sort(width.v[!duplicated(width.v)], decreasing=T)
  
  print(paste0(chr, ": Widths are:" ), quote=F)
  print(width.v, quote=F)
  
  # Identify persistent unmappable areas
  temp <- c(paste0("x.", ct.v), paste0("y.", ct.v))
  persist.TF <- by(data=paste0(x$axis, ".", x$ct),
                   INDICES=sort.id, FUN=function(x){
                     identical(sort(x), sort(temp))
                   })
  
  # Area ids unmappable across all tissues
  unmap.bin <- names(persist.TF)[as.vector(persist.TF)]
  
  # Largest unmappable area
  unmapbig.bin <- names(width.v)[width.v==max(width.v)]
  
  if( length(unmapbig.bin)!=1 ){
    stop(paste0(chr, ": Multiple largest unmappable area."))
  }
  
  if(!unmapbig.bin%in%unmap.bin){
    stop(paste0(chr, ": Largest unmappable not unmappable across all tissues"))
  } 
  
  # Largest unmappable area
  eval(parse(text=paste0(
    "unmapbig.bin <- c(", unmapbig.bin, ")"
  )))
  unmapbig[[chr]] <- c(chr, min(unmapbig.bin), max(unmapbig.bin))
  
  # Rest of persistent unmappable bins
  unmap.bin <- strsplit(x=unmap.bin, split=":")
  unmap.bin <- do.call(rbind, unmap.bin)
  unmap[[chr]] <- c(chr, paste(unmap.bin[,1], collapse=";"), paste(unmap.bin[,2], collapse=";"))
  
  rm(unmap.bin, unmapbig.bin, sort.id, x, data.name, width.v, temp, persist.TF)
  
}

unmap <- do.call("rbind.data.frame", c(unmap, stringsAsFactors=FALSE))
colnames(unmap) <- c("chr", "start", "end")
write.csv(unmap, file=paste0(out.dir, "/", out.name, "_unmappable.csv"), 
          row.names=FALSE)
rm(unmap)

unmapbig <- do.call("rbind.data.frame", c(unmapbig, stringsAsFactors=FALSE))
colnames(unmapbig) <- c("chr", "start", "end")
write.csv(unmapbig, file=paste0(out.dir, "/", out.name, "_unmappablebig.csv"), 
          row.names=FALSE)

# rm(list=ls()); gc()