################################################################################
# Plot mean Cp (Cp=0:21) vs. mean Cp (Cp=1:21)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/2_Expl_contacts"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
meanCp.id = "Cp1To21"
meanCp.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/out_binWeightedMeanCp/Cp1To21")
out.dir  = paste0(wk.dir, "/out_meanCp0vsMeanCp1/", meanCp.id)
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X")) 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/lmEqn_string.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( !dir.exists(out.dir) ){
  print("Creating output directory...", quote=FALSE)
  dir.create(path=out.dir)
}

df <- data.frame()
for(chr in chr.v){
  
  load(file=paste0(meanCp.dir, "/", gcb, "_", chr, "_weightedMeanCp.RData"))
  
  df <- rbind(df, BINWMEANCP.DF[,colnames(BINWMEANCP.DF)!="bin"])
  rm(BINWMEANCP.DF); gc()
  
  print(paste0(chr, " done!"), quote=FALSE)
  
}
colnames(df)[colnames(df)==c("wmeanCp0")] <- "meanCp0"
colnames(df)[colnames(df)==c("wmeanCp")] <- "meanCp1"

incl.TF <- is.finite(df$meanCp0) & is.finite(df$meanCp1)
df <- df[incl.TF,]

out.id <- paste0(gcb, "_meanCpID", meanCp.id)
p.title <- paste0(out.id, "_", length(incl.TF), "origbins_", sum(incl.TF), 
                  "binsWithMeanCp0ANDCp1")
rm(incl.TF)

# meanCp0 vs. meanCp1

y <- "meanCp0"

lm.str <- lmEqn_string(x="meanCp1", y=y, data=df)

p <- ggplot(data=df, aes_string(x="meanCp1", y=y)) +
  geom_point(alpha=0.05) +
  labs(title=p.title) + 
  annotate(geom="text", x=max(df$meanCp1), y=min(df[[y]]),
           size=8, parse=TRUE, hjust=1, vjust=0,
           label=lm.str) +
  bgr2

ggsave(filename=paste0(out.dir, "/", out.id, "_", y, "VSmeanCp1_scat.pdf"), 
       plot=p, height=10, width=10)

p <- ggplot(data=df, aes_string(x="meanCp1", y=y)) +
  geom_hex() +
  labs(title=p.title) + 
  annotate(geom="text", x=max(df$meanCp1), y=min(df[[y]]),
           size=8, parse=TRUE, hjust=1, vjust=0,
           label=lm.str) +
  bgr2

ggsave(filename=paste0(out.dir, "/", out.id, "_", y, "VSmeanCp1_hex.pdf"), 
       plot=p, height=10, width=10)

df$LRijDIVposLRij <- df$LRij/df$posLRij

# meanCp0/1 vs. LRijDIVposLRij

for( y in c("meanCp0", "meanCp1") ){
  
  lm.str <- lmEqn_string(x="LRijDIVposLRij", y=y, data=df)
  
  p <- ggplot(data=df, aes_string(x="meanCp1", y=y)) +
    geom_point(alpha=0.05) +
    labs(title=p.title) + 
    annotate(geom="text", x=max(df$LRijDIVposLRij), y=min(df[[y]]),
             size=8, parse=TRUE, hjust=1, vjust=0,
             label=lm.str) +
    bgr2
  
  ggsave(filename=paste0(out.dir, "/", out.id, "_", y, "VSLRijDIVposLRij_scat.pdf"), 
         plot=p, height=10, width=10)
  
  p <- ggplot(data=df, aes_string(x="LRijDIVposLRij", y=y)) +
    geom_hex() +
    labs(title=p.title) + 
    annotate(geom="text", x=max(df$LRijDIVposLRij), y=min(df[[y]]),
             size=8, parse=TRUE, hjust=1, vjust=0,
             label=lm.str) +
    bgr2
  
  ggsave(filename=paste0(out.dir, "/", out.id, "_", y, "VSLRijDIVposLRij_hex.pdf"), 
         plot=p, height=10, width=10)
  
}

# rm(list=ls()); gc()