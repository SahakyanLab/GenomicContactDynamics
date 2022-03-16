################################################################################
# Convert depleted, priority foi text files into RData with feature names
################################################################################
wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/7_FeaturePermutation/out_table/feat_repeats"
depleted.file = paste0(wk.dir, "/foifile_depleted_nperm10000_seed834_mxmskfr0_CptopCP3_pvalcutoff0.05_numOlapA")
priority.file = paste0(wk.dir, "/foifile_priority_nperm10000_seed834_mxmskfr0_CptopCP3_pvalcutoff0.05_numOlapA")
out.name = "foifile_depleted_priority_nperm10000_seed834_mxmskfr0_CptopCP3_pvalcutoff0.05_numOlapA"

FOI <- list()
for( x in c("depleted", "priority") ){
  
  eval(parse(text=paste0(
    'foi.v <- readLines(con=paste0(', x, '.file))'
  )))
  
  
  foi.v <- strsplit(x=foi.v, split="_foi_|_desc_")
  foi.v <- unlist(lapply(X=foi.v, FUN=function(x) x[2]))
  
  FOI[[x]] <- foi.v
  
}

save(x=FOI, file=paste0(wk.dir, "/", out.name, ".RData"))

# rm(list=ls()); gc()