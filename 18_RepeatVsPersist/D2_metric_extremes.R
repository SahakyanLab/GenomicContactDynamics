################################################################################
# Per subfamily, determine proportion of contacts with at least 1 shared number
# of site considering only contacts with at least 2 sites of the subfamily so as 
# to level the ground between subfamilies with different copy numbers (MINREPCOUNTS). 
# Using Cps.forCalc, specify which Cps to consider for calculating proportion.
# If Cps.forCalc pertains to more than 1 Cp value, the proportion in final data
# would be the mean of the proportions from the Cp values in Cps.forCalc. 
# Append subfamily groupings for further characterisations. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")

metric = "minrep"
src.dir = paste0(wk.dir, "/z_ignore_git/out_minRepCounts/subfamALL_", metric, "_atleast2sumrep")
out.dir = paste0(wk.dir, "/out_combine") #paste0(wk.dir, "/out_metric_extremes")

repeatmx.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/17_RepeatAge/z_ignore_git/out_addToSummary")
repeat.file = paste0(repeatmx.dir, "/hg19repeats_repName.RData")
#element.file = paste0(wk.dir, "/Repeat_rankingbyAge/plot_GiorPubl372rankrepFamilies.csv")
### OTHER SETTINGS #############################################################
src.nme = paste0("chrALL_min2Mb_subfamALL_", metric, "Counts")
Cps.forCalc = 19:21 # Cp values to consider when calculating median non-zero count fraction
out.id = paste0("metric_extremes_Cp", Cps.forCalc[[1]], "To", tail(Cps.forCalc, n=1))
GROUP.CLASS <- list(
  not.transposon = c("Low_complexity", "RNA", "rRNA", "Satellite", "scRNA", 
                     "Simple_repeat", "snRNA", "srpRNA", "tRNA"),
  DNA.transposon = c("DNA", "DNA?"), 
  retro.transposon = c("LINE", "SINE", "LTR", "RC", "LINE?", "SINE?", "LTR?"),
  not.classified = c("Other", "Unknown", "Unknown?")
)
group.cols = setNames(object=c("#7E6148FF", "#3C5488FF", "#00A087FF", "black"), 
                      nm=names(GROUP.CLASS))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
#library(ggsci)
#library(ggplot2)
#source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(paste0(src.dir, "/", src.nme, ".RData"))
elements <- names(MINREPCOUNTS)

# Deal with NULLS, these are for Cp values without contacts with at least 2 sites
# of the element; convert nulls to c(`0`=NA)
for(elm in elements){
  
  MINREPCOUNTS[[elm]] <- sapply(names(MINREPCOUNTS[[elm]]), simplify=F, FUN=function(cp.nme){
    counts <- MINREPCOUNTS[[elm]][[cp.nme]]
    if(is.null(counts)){
      return(c(`0`=NA)) 
    } else {
      return(counts)
    }
  })
  
}

if( unique(lengths(MINREPCOUNTS)) != 21 ){
  stop("Missing Cps.")
}

# Sum of contact counts combining Cps.forCalc

SUM.COUNTS <- MINREPCOUNTS
for(elm in elements){
  SUM.COUNTS[[elm]] <- sum(unlist( MINREPCOUNTS[[elm]][ Cps.forCalc ] ), na.rm=T)
  #SUM.COUNTS[[elm]] <- lapply(MINREPCOUNTS[[elm]], FUN=sum, na.rm=T)
}

if( any( !is.finite(unlist(SUM.COUNTS)) ) ){
  stop("Non-finite sum.")
}

# Contact counts with 0 value combining Cps.forCalc

NONZERO.COUNTS <- MINREPCOUNTS
for(elm in elements){
  
  nonzero.counts <- lapply(MINREPCOUNTS[[elm]][ Cps.forCalc ], FUN=function(cp.counts){
    cp.counts <- cp.counts[names(cp.counts) != "0"]
    return(cp.counts)
  })
  NONZERO.COUNTS[[elm]] <- sum(unlist( nonzero.counts ), na.rm=T)
  
  #NONZERO.COUNTS[[elm]] <- lapply(MINREPCOUNTS[[elm]][], FUN=function(counts){
  #  counts <- counts[names(counts) != "0"]
  #  return(sum(counts))
  #})

}

# Fraction of contacts with 0s 

NONZERO.COUNTS.FR <- MINREPCOUNTS
for(elm in elements){
  
  val <- NONZERO.COUNTS[[elm]]
  ref <- SUM.COUNTS[[elm]]
  
  NONZERO.COUNTS.FR[[elm]] <- val / ref
  
  #val <- unlist(NONZERO.COUNTS[[elm]])
  #ref <- unlist(SUM.COUNTS[[elm]])
  
  #if( identical(names(val), names(ref)) ){
  #  NONZERO.COUNTS.FR[[elm]] <- as.list(val / ref)
  #} else {
  #  stop(paste0(elm, ": Cp names not matching."))
  #}
  
}

## Min fraction of contacts with 0s across Cps per element

#NONZERO.COUNTS.FR.PLOT <- MINREPCOUNTS
#for(elm in elements){
#  NONZERO.COUNTS.FR.PLOT[[elm]] <- mean(unlist( NONZERO.COUNTS.FR[[elm]][as.character(Cps.forCalc)] ))
#}

## PLOT data

df <- stack(unlist(NONZERO.COUNTS.FR))
df$ind <- as.character(df$ind)
setnames(df, old="ind", new="repName")

# Add subfamily info
load(file=repeat.file)
rownames(REPEAT.MX) <- NULL
REPEAT.MX <- REPEAT.MX[,c("repName", "repClass", "repFamily")]

df <- merge(x=df, y=REPEAT.MX, by="repName", all=T)

grp.class.df <- stack(GROUP.CLASS)
grp.class.df$ind <- as.character(grp.class.df$ind)
setnames(grp.class.df , old=c("values", "ind"), new=c("repClass", "plot.group"))

df <- merge(x=df, y=grp.class.df, by="repClass", all.x=T)
rm(grp.class.df)

#

df <- df[order(df$repName, decreasing=F),]
rownames(df) <- NULL
save(df, file=paste0(out.dir, "/", metric, "_", out.id, ".RData"))

# ## PLOTS
# 
# df <- df[order(df$values, decreasing=T),]
# df$plot.group <- factor(df$plot.group, levels=names(GROUP.CLASS))
# df$repName <- factor(df$repName, levels=unique(df$repName))
# 
# p <- ggplot(data=df, aes(x=repName, y=values)) +
#   #geom_col(aes(col=plot.group, fill=plot.group)) +
#   geom_segment(aes(x=repName, xend=repName, y=0, yend=values, col=plot.group), alpha=0.2) +
#   geom_point(aes(col=plot.group), shape=1, size=2) + 
#   #geom_hline(yintercept=0.5, linetype="dashed", color="black", size=0.5) + 
#   scale_color_manual(values=group.cols) + 
#   bgr2 + 
#   theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
#         axis.line.x=element_line(color="black"), axis.line.y=element_blank(),
#         panel.background=element_blank(), legend.position="none") +
#   facet_grid(.~plot.group) 
# 
# ggsave(filename=paste0(out.dir, "/", src.nme, "_", out.id, ".png"),
#        width=10*300, height=10*300, unit="px", plot=p)

# rm(list=ls()); gc()