################################################################################
# Identify potential unmappable areas, areas that don't form contacts in all 
# tissues.
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
data.dir = paste0(wk.dir, "/out_contactMxGap")
out.dir = paste0(wk.dir, "/out_unmappable")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))
metric = "Cs.norm"
out.id = "whole"
ct.v = sort(c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
              "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC"))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
ct.v.len <- length(ct.v)
out.name <- paste(gcb, out.id, metric, sep="_")

unmap <- list()
for(chr in chr.v){
  
  data.name <- paste(gcb, chr, out.id, metric, sep="_")
  x <- read.csv(file=paste0(data.dir, "/", data.name, "_contactMxGap.csv"), 
                row.names="X")
  if( any(duplicated(x)) ){
    stop(paste0(chr, " duplicated rows."))
  }
  
  sort.id <- paste0(x$start, ":", x$end)
  count.id <- table(sort.id)
  unmap.bin <- names(count.id)[count.id==(ct.v.len*2)]
  unmap.bin <- paste(unmap.bin, collapse=",")
  unmap.bin <- eval(parse(text=paste0(
    "unmap.bin <- c(", unmap.bin, ")"
  )))
  unmap.bin <- sort(unique(unmap.bin))
  unmap[[chr]] <- c(chr, paste(unmap.bin, collapse=";"))
  
  rm(unmap.bin, sort.id, count.id, x, data.name)
}

unmap <- do.call("rbind.data.frame", c(unmap, stringsAsFactors=FALSE))
colnames(unmap) <- c("chr", "unmapbin")
write.csv(unmap, file=paste0(out.dir, "/", out.name, "_unmappable.csv"), 
          row.names=FALSE)

# rm(list=ls()); gc()