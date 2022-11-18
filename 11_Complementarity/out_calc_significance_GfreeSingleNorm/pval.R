################################################################################
# Get which Cp pairs are not significant i.e. p-value >= alpha
################################################################################
library(reshape2)
alpha = 0.05

# In df, df$Var1 > df$Var2
Cp.select1 = 17:21
Cp.select2 = 1:5

pw.dir = "~/SahakyanLab/GenomicContactDynamics/11_Complementarity/out_calc_significance_GfreeSingleNorm/filter_gap_jMINUSiMINUS1equals50To100bins" #/filter_gap_jMINUSiMINUS1equals50bins"
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

#

is.Cp.select <- df$Var1 %in% Cp.select1 & df$Var2 %in% Cp.select2
df <- df[df$value >= alpha & is.Cp.select,]

print(df)

# rm(list=ls()); gc()