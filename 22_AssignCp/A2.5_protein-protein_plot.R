################################################################################
# Plot protein-protein metrics vs. Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Avoid left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/X1_")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/X1_")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
data.file = paste0(wk.dir, "/out_protein-protein/min2Mb_protein-protein_maxCp_full.csv")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/compareTwoDist.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
df <- data.table::fread(file=data.file, stringsAsFactors=F, header=T, data.table=F, check.names=T)
df <- df[!is.na(df$Cp),]
df$Cp <- factor(as.character(df$Cp), levels=as.character(1:21))

#> table(df$Cp)

#1     2     3     4     5     6     7     8     9    10    11    12    13 
#43532 41186 35528 28734 22219 17641 13438 10661  7917  6183  4766  3577  2627 
#14    15    16    17    18    19    20    21 
#2084  1463  1030   662   420   228   106    24 

vals <- setdiff(colnames(df), c("X4.Gene1", "X5.Gene2", "chrom", "Cp", "X2.Gene1", "X3.Gene2",
                                "X15.LLR_GIN", "X24.LLR_CIN", # These columns, all NAs
                                "X39.Max_category"))

out.name <- gsub(data.file, pattern=".csv", replacement="_plot")

# Density plot: Cp >= 3 vs. Cp <= 19 distributions, perform significance

df$group <- NA
df$group[as.character(df$Cp) %in% c("1", "2", "3")] <- "Cp1To3"
df$group[as.character(df$Cp) %in% c("19", "20", "21")] <- "Cp19To21"
df$group <- factor(df$group, levels=c("Cp1To3", "Cp19To21"))
cols <- c("#5E4FA2FF", "#9E0142FF")

is.with.group <- !is.na(as.character(df$group))
is.Cp1To3 <- as.character(df$group) == "Cp1To3"
is.Cp19To21 <- as.character(df$group) == "Cp19To21"
  
p.lst <- list()
for(val in vals){
  
  is.valid <- is.nonNAval & is.with.group
 
  TEST <- compareTwoDist(x=df[[val]][is.valid & is.Cp1To3],
                         y=df[[val]][is.valid & is.Cp19To21])
  
  p.lst[[val]] <- ggplot(data=df[is.valid,], aes_string(x=val)) +
    geom_density(aes(fill=group, col=group), size=1.5, alpha=0.7) + 
    scale_color_manual(values=cols) + 
    scale_fill_manual(values=cols) +
    labs(title=paste0(TEST$test.id, "\n x=Cp1To3, y=Cp19To21"), y=NULL, x=val) +
    bgr2 +
    theme(axis.title.x=element_text(size=10),
          plot.title=element_text(size=4))
  
}

p.arr <- ggarrange(plotlist=p.lst, nrow=6, ncol=6, common.legend=T)
ggexport(p.arr, height=30, width=30, filename=paste0(out.name, "_Cp1To3VsCp19To21.pdf"))

# Box plot

png(filename=paste0(out.name, ".png"), height=100*20, width=100*40, res=100)
par(mfrow=c(6,6))

for(val in vals){
  
  setnames(x=df, old=val, new="val")
  boxplot(formula=val~Cp, data=df, outline=F, xlab="Cp", ylab=val, 
          main="", col="#FDC776", cex.main=0.8)
  df$val <- NULL
  
  print(val, quote=F)
  
}

dev.off()

# rm(list=ls()); gc()