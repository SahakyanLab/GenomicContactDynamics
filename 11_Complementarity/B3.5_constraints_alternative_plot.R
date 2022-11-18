################################################################################
# Alternative plot for complementarity values
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

options(warnPartialMatchDollar = T)
options(warn = 1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/11_Complementarity")
lib = paste0(home.dir, "/DPhil/lib")
constraints.id = "GfreeSingleNorm" #"hg19_rm_GfreeSingleNorm"
compl.dir = paste0(wk.dir, "/out_constraints_", constraints.id, "/merged_final")
out.dir = paste0(wk.dir, "/out_constraints_alternative_plot")
### OTHER SETTINGS #############################################################
chr = "chrALL" 
gcb = "min2Mb"
type = "kmer"
affix = ""
Cps = 1:21
ylim.val = list(CIIkmer=c(-2.5,0), CIIalign=NULL, 
                Gfreekmer=c(-1.2,-0.4), sdDifferencekmer=c(0,0.00030))

out.id = constraints.id # gap_effect

DATA.PATH = list(paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"),
                 paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"),
                 paste0(compl.dir, "/", chr, "_", type, "_", gcb, affix, ".RData"))
GAP.RNG = list(NULL, c(50,100), c(50,50)) # j - i - 1, closed range, set to NULL if not filtering
lty.val = c("solid", "dashed", "dotted")
shp.val = c(15,1,2)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
LENS <- c(length(DATA.PATH), length(GAP.RNG), length(lty.val), length(shp.val))
if( any(LENS != 3) ){
  rm(GAP.RNG, DATA.PATH)
  stop("Check lengths of objects, whose lengths should be equal.")
}

if( !is.null(GAP.RNG) ){
  GAP.RNG.ID <- lapply(X=GAP.RNG, FUN=function(gap.rng){
    if(is.null(gap.rng)){
      return("NULLgap")
    } else {
      return( paste0("gap", paste(gap.rng, collapse="To"), "bins") )
    }
  })
}

#

df <- sapply(X=1:length(DATA.PATH), simplify=F, FUN=function(ind){
  
  load(DATA.PATH[[ind]])
  gap.rng <- GAP.RNG[[ind]]
  if( !is.null(gap.rng) ){
    print(paste0(GAP.RNG.ID[[ind]], ": Filtering contacts based on this closed gap range."))
    gaps <- CII.MX[,"j"] - CII.MX[, "i"] - 1
    CII.MX <- CII.MX[ gaps >= gap.rng[[1]] & gaps <= gap.rng[[2]], ]
    rm(gaps)
  } else {
    print(paste0(GAP.RNG.ID[[ind]], ": No filtering of contacts based on gap."))
  }
  CII.MX <- data.frame(CII.MX, dta.id=GAP.RNG.ID[[ind]], check.names=F)
  return(CII.MX)
  
})

# Plot

df <- do.call("rbind", df)
setnames(x=df, old="C||", new="CII")

# Remove contacts with no Cp
df <- df[!is.na(df$Cp),]

df$Cp <- factor(as.character(df$Cp), levels=as.character(Cps))
df$dta.id <- factor(as.character(df$dta.id), levels=unlist(GAP.RNG.ID))

compl.types <- setdiff(colnames(df), c("i", "j", "Cp", "dta.id"))
out.name <- paste0(out.id, "_", chr, "_", type, "_", gcb, affix)
dta.ids <- levels(df$dta.id)
gaps.noNACp <- df$j - df$i - 1 

for(c.type in compl.types){
  
  # Calculate actual range of gap values and add to title
  
  is.nonNA.c.type <- !is.na(df[[c.type]])
  
  ACT.GAP.RNG.ID <- sapply(X=dta.ids, simplify=T, FUN=function(dta.id){
    is.dta <- df$dta.id == dta.id
    gap.rng.noNA <- range( gaps.noNACp[is.nonNA.c.type & is.dta] ) 
    return( paste0("gap", paste(gap.rng.noNA, collapse="To"), "bins") )
  })
 
  plot.title <- paste0(out.name, "\nactualgapAfterNAxyvaldropped_", paste(ACT.GAP.RNG.ID, collapse="_"))
  
  p <- ggplot(data=df[is.nonNA.c.type,], aes_string(x="Cp", y=c.type)) +
    stat_boxplot(col="gray50", geom="errorbar", aes(lty=dta.id), width=0, position=position_dodge(0.7)) + # To add only whiskers of default boxplot
    stat_summary(col="#fc9a08", fun="median", aes(shape=dta.id), position=position_dodge(0.7)) +
    scale_y_continuous(limits=ylim.val[[paste0(c.type, type)]]) + 
    scale_linetype_manual(values=lty.val) + 
    scale_shape_manual(values=shp.val) + 
    labs(title=paste0(plot.title, "\ngapequalsjMINUSiMINUS1_pointsAtMedian_errorbarSameAsDefaultBoxplotWhiskers")) +
    bgr2 +
    theme(aspect.ratio=0.8,
          legend.position="bottom",
          plot.title=element_text(size=5))
  
  ggsave(filename=paste0(out.dir, "/", out.name, "_", c.type, ".pdf"), 
                         width=10, height=8, plot=p)
  
  rm(is.nonNA.c.type, ACT.GAP.RNG.ID, plot.title, p)
  gc()
  
  print(paste0(c.type, " done!"), quote=F)
  
}

# rm(list=ls()); gc()






