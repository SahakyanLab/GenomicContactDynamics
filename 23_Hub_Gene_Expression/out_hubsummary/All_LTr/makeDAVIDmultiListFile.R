################################################################################
# Make DAVID-compatible multi-list file.

wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer/out_hubsummary/All_LTr"
out.dir = wk.dir

fle.v <- list.files(path=wk.dir, pattern="_genes.txt", full.names=F)
if( length(fle.v)!= 38){
  stop("Wrong number of files.")
}

MLIST.DF <- list()
for(fle in fle.v){
  
  MLIST.DF[[fle]] <- readLines(con=paste0(wk.dir, "/", fle))
   
}

ngenes <- unlist(lapply(X=MLIST.DF, FUN=length))
ngenes.max <- max(ngenes)

MLIST.DF <- lapply(X=MLIST.DF, FUN=function(hubgenes){
  
  add.len <- ngenes.max-length(hubgenes)
  c(hubgenes, rep(x="", times=add.len))
  
})

MLIST.DF <- do.call("cbind", MLIST.DF)

write.table(x=MLIST.DF, file=paste0(out.dir, "/min2Mb_All_topCP3_gapBin50_multilist.txt"),
            sep="\t", row.names=F)

# rm(list=ls()); gc()
