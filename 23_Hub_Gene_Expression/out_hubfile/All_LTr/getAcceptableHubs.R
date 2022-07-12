################################################################################
wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer/out_hubfile"
flenme.ref.dir = paste0(wk.dir, "/All_AllTr/acceptable_n")

flesrc.dir = paste0(wk.dir, "/All_LTr")
out.dir = paste0(flesrc.dir, "/acceptable_n")

flenme.v <- list.files(path=flenme.ref.dir, full.names=F)
flenme.v <- flenme.v[flenme.v != "min2Mb_chr1_All_topCP3_gapBin50_hub4.csv"]
rm(flenme.ref.dir)

for(fle in flenme.v){
  
  file.copy(from=paste0(flesrc.dir, "/", fle), to=out.dir)
  print(fle, quote=F)
  
}

# rm(list=ls()); gc()