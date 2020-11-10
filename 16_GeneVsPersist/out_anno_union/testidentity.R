dir.old <- "/Users/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist/out_anno_union"
dir.new <- "/Users/ltamon/Desktop/out_anno_union"

f.old <- list.files(dir.old)
f.new <- list.files(dir.new)

for(f in f.old){
  old <- readLines(paste0(dir.old, "/", f))
  new <- readLines(paste0(dir.new, "/", f))
  if( !identical(old, new) ){
    stop(f)
  }
  print(identical(old, new))
}
