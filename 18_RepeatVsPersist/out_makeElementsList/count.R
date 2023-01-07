################################################################################
# Sum count of elements per chr to check if sum matches with MINELM.MX SOURCE
# generated.
################################################################################
wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist/out_makeElementsList/subfamALL/"
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X"))

el.len <- 0
for(chr in chr.v){
  
  el.v <- readLines(con=paste0(wk.dir, "/", gcb, "_", chr, "_elements.txt"))
  el.len <- el.len + length(el.v)
  
  print(chr, quote=FALSE)
  
}
print(paste0(el.len, " in total."), quote=FALSE)
################################################################################
# rm(list=ls()); gc()