################################################################################
# Get which Cp pairs are not significant i.e. p-value >= alpha
################################################################################
library(reshape2)
alpha = 0.01

pw.dir = "~/SahakyanLab/GenomicContactDynamics/11_Complementarity/out_calc_significance_GfreeSingleNorm"
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