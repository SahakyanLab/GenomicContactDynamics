################################################################################
# Per fature-based contact type, plot fraction of contacts across Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=TRUE)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon/"
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib      = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir   = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/7_FeaturePermutation")

featvscp.dir = paste0(wk.dir, "/out_foiVsij")
out.dir = paste0(wk.dir, "/out_frijPrCp")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
olap.type = "any_querybin" # "within_querybin" | "within_queryfoi"
foi.num.v = 194
# To distinguish generic names e.g. none-none, append foi.num
foi.pair.v = c("hg19_TOP2B_TF_MCF7-hg19_TOP2B_TF_MCF7.194.",
               "none-hg19_TOP2B_TF_MCF7.194.") #,
#foi.pair.v = c("none-hg19_TOP2B_TF_MCF7.194.")#,
#               "ESC_TOP2B_TF_MCF7-ESC_TOP2B_TF_MCF7.195.",
#               "none-ESC_TOP2B_TF_MCF7.195.",
#               "FC_TOP2B_TF_MCF7-FC_TOP2B_TF_MCF7.196.",
#               "none-FC_TOP2B_TF_MCF7.196.",
#               "LC_TOP2B_TF_MCF7-LC_TOP2B_TF_MCF7.197.",
#              "none-LC_TOP2B_TF_MCF7.197.")
combinePairs = T # If T, add counts 
col.v = "gray40" #c("#00A087FF", "#F39B7FFF") # Order should correspond to foi.pair.v
#col.v = "#F39B7FFF"
out.id = "hg19_TOP2B-TOP2B_TF_MCF7_hg19_none-TOP2B_TF_MCF7"
#out.id = "hg19_TOP2B-TOP2B_TF_MCF7"
#out.id = "hg19_none-TOP2B_TF_MCF7"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
foi.num.v <- as.character(foi.num.v)

FEATVSCP <- NULL
for(foi.num in foi.num.v){
  
  fle <- list.files(path=featvscp.dir, pattern=paste0("_", foi.num, "_", olap.type, 
                                                      "_foiVsij.RData"), full.names=F)
  if(length(fle)!=1){
    stop(paste0(foi.num, ": Multiple .RData"))
  } else {
    
    load(file=paste0(featvscp.dir, "/", fle))
    
    allijPCp <- FEATVSCP.MX[, "All"]
    FEATVSCP.MX <- FEATVSCP.MX[, dimnames(FEATVSCP.MX)[[2]]!="All"]
    
    ## Convert counts of contacts to fraction per Cp
    #FEATVSCP.MX <- FEATVSCP.MX/allijPCp
    
    # Removed this checkpoint because of precision problems with float
    #if( any(rowSums(FEATVSCP.MX)!=1) ){
    #  stop(paste0(foi.num, ": Fractions do not sum to 1."))
    #} else {
      
      # To distinguish generic names e.g. none-none, append foi.num
      dimnames(FEATVSCP.MX)[[2]] <- paste0(dimnames(FEATVSCP.MX)[[2]], ".", foi.num, ".")
      
      FEATVSCP <- cbind(FEATVSCP, FEATVSCP.MX)
      
    #}
    
  }
  
  rm(FEATVSCP.MX)
  print(paste0(foi.num, " done!"), quote=F)
  
}

foi.pair1.v <- intersect(foi.pair.v, dimnames(FEATVSCP)[[2]])
if( !identical(foi.pair.v, foi.pair1.v) ){
  stop("Element/s in foi.pair.v not in FEATVSCP.")
} else {
  FEATVSCP <- FEATVSCP[,foi.pair.v, drop=F]
}
rm(foi.pair1.v)

ncol.tmp <- ncol(FEATVSCP)
if(combinePairs){
  
  FEATVSCP <- as.matrix(rowSums(x=FEATVSCP, na.rm=F), ncol=ncol.tmp)
  foi.pair.v <- paste(foi.pair.v, collapse="_")
  dimnames(FEATVSCP)[[2]] <- foi.pair.v
  
}

FEATVSCP <- FEATVSCP/allijPCp
  
# Plot

Cp.levels <- as.character(dimnames(FEATVSCP)[[1]])
Cp.levels <- as.character(sort(as.numeric(Cp.levels), decreasing=F))
FEATVSCP <- as.data.frame(FEATVSCP)
FEATVSCP$Cp <- factor(x=as.character(dimnames(FEATVSCP)[[1]]),
                      levels=Cp.levels)
rm(Cp.levels)

df <- reshape2::melt(data=FEATVSCP, id="Cp")
df$variable <- as.character(df$variable)
df$variable <- factor(x=df$variable, 
                      levels=intersect(foi.pair.v, unique(df$variable)))

out.name <- paste0(gcb, "_", olap.type, "_", out.id, "_combinePairs", combinePairs)
p <- ggplot(data=df, aes(x=Cp, y=value)) + 
  geom_point(aes(colour=variable), size=5) +
  #ggsci::scale_colour_npg() + 
  scale_colour_manual(values=col.v) + 
  labs(y="Contact fraction", title=out.name) + 
  bgr2 +
  theme(legend.text=element_text(size=1),
        aspect.ratio=0.40)

ggsave(filename=paste0(out.dir, "/", out.name, "_frijVsCp_scatter.pdf"),
       width=10, height=5, plot=p)


# Fold change
df$value <- log2(df$value/df$value[df$Cp==1])
p <- ggplot(data=df, aes(x=Cp, y=value)) + 
  geom_point(aes(colour=variable), size=5) +
  #ggsci::scale_colour_npg() + 
  scale_colour_manual(values=col.v) + 
  labs(y="Fold change", title=out.name) + 
  bgr2 +
  theme(legend.text=element_text(size=1),
        aspect.ratio=0.40)

ggsave(filename=paste0(out.dir, "/", out.name, "_foldchange_frijVsCp_scatter.pdf"),
       width=10, height=5, plot=p)

# rm(list=ls()); gc()


