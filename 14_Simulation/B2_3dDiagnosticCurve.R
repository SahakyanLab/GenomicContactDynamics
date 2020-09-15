################################################################################
# ROC 
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
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
csv.dir = paste0(wk.dir, "/out_compare")
out.dir = paste0(wk.dir, "/out_curve")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr1"
ct = "FC"
metric.v = c(subj="SIM.4.2.kmer.5", ref="Cs.norm")
out.id = "only_gap50To1900_inclu1000To3200" 

# Closed ranges where most contacts belong; to colour points differently
subj.range = c(0,0.0015)
ref.range = c(0,0.4)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(plotly)
source(paste0(wk.dir, "/lib/makeInteractiveScatter.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
id <- paste(paste(names(metric.v), metric.v, sep="_"), collapse="_")
out.name <- paste(gcb, chr, ct, out.id, id, sep="_"); rm(id)
print(paste0(out.name, "..."), quote=FALSE)

COMPIJMX <- read.csv(file=paste0(csv.dir, "/", out.name, ".csv"),
                     header=TRUE, stringsAsFactors=FALSE)

TPR <- COMPIJMX$TP/(COMPIJMX$TP+COMPIJMX$FN)
FPR <- COMPIJMX$FP/(COMPIJMX$FP+COMPIJMX$TN)
prec <- COMPIJMX$TP/(COMPIJMX$TP+COMPIJMX$FP)

#c.off <- COMPIJMX$c.offsubj<0 | COMPIJMX$c.offref<0
#incl.TF <- is.finite(TPR) & is.finite(FPR) & is.finite(prec) #& c.off
#incl.TF <- rep(TRUE, nrow(COMPIJMX))

out.range <- (COMPIJMX$c.offsubj<subj.range[1] | COMPIJMX$c.offsubj>subj.range[2]) |
  (COMPIJMX$c.offref<ref.range[1] | COMPIJMX$c.offref>ref.range[2])
group <- out.range
group[group==TRUE] <- "out-range"
group[group==FALSE] <- "in-range"

df <- cbind.data.frame(c.offsubj=COMPIJMX$c.offsubj, c.offref=COMPIJMX$c.offref,
                       TPR=TPR, FPR=FPR, prec=prec, cutoff=rownames(COMPIJMX),
                       cutoff.lab=paste0("Cut-off ind:", rownames(COMPIJMX), 
                                         " subj:", COMPIJMX$c.offsubj, 
                                         " ref:", COMPIJMX$c.offref),
                       group=group); rm(group)
col.v <- c(`in-range`="darkblue", `out-range`="gray")[levels(df$group)]
#df <- df[incl.TF,]

# Precisio-recall(TPR) plot
makeInteractiveScatter(
  out.dir=out.dir,
  out.name=paste0(out.name, "_precall"),
  df=df,
  x="TPR",
  y="prec",
  z="cutoff",
  grp="group",
  xlim.v=NULL,
  ylim.v=NULL,
  zlim.v=NULL,
  point.lab="cutoff.lab",
  col.v=col.v,
  main=paste0(out.name, "\n x=TPR y=prec z=cut-off ind")
)

# ROC
makeInteractiveScatter(
  out.dir=out.dir,
  out.name=paste0(out.name, "_roc"),
  df=df,
  x="FPR",
  y="TPR",
  z="cutoff",
  grp="group",
  xlim.v=NULL,
  ylim.v=NULL,
  zlim.v=NULL,
  point.lab="cutoff.lab",
  col.v=col.v,
  main=paste0(out.name, "\n x=FPR y=TPR z=cut-off ind")
)
# Combined
makeInteractiveScatter(
  out.dir=out.dir,
  out.name=paste0(out.name, "_combined"),
  df=df,
  x="prec",
  y="FPR",
  z="TPR",
  grp="group",
  xlim.v=NULL,
  ylim.v=NULL,
  zlim.v=NULL,
  point.lab="cutoff.lab",
  col.v=col.v,
  main=paste0(out.name, "\n x=prec y=FPR z=TPR")
)

# rm(list=ls()); gc()
