################################################################################
# Get which Cp pairs are not significant i.e. p-value >= alpha
################################################################################
library(reshape2)
alpha = 0.05

pw.dir = "~/SahakyanLab/GenomicContactDynamics/20_ReplicationTiming/out_contact_Cp_plot/chrALL"
pw.files <- list.files(pw.dir, pattern="pairwisedifftest.RData")

df <- list()
for(pw in pw.files){
  
  load(paste0(pw.dir, "/", pw))
  TEST$meanmed <- TEST$alt <- NULL
  
  TEST$pt <- melt(TEST$pt$p.value)
  TEST$pmw <- melt(TEST$pmw$p.value)
  
  df[[pw]] <- do.call("rbind", TEST)

  print(pw, quote=F)
  
  rm(TEST, pw)
  
}

df <- do.call("rbind", df)
df <- df[!is.na(df$value),]
df <- df[df$value >= alpha,]

print(df)

# rm(list=ls()); gc()

rw.nme <- gsub(pattern="min2Mb_allChr_ijratesVsCp_min.countPerBin3_ij_recomRates_2011_01_phaseII_B37_Myers_Cp0To21compare_pairwisedifftest.RData",
               replacement="", x=rownames(df))
rownames(df) <- rw.nme
