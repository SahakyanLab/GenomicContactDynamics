################################################################################
# Calculate confusion matrix-derived metrics
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
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
csv.dir = paste0(wk.dir, "/out_compare/csv_chr1_whole_maskMidSquare_gap50up_refCp")
out.dir = paste0(wk.dir, "/out_metricVolume")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr = "chr1" 
ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
         "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")

# Subj should have the same scale of cut-off values
Metric = list(subj=c("SIM.4.2.kmer.5", "SIM.3.2.kmer.5"), 
              # One ref only
              ref="Cp")
src.id = "whole_maskMidSquare_gap50up"
out.id = "SIMVsCp_partRange"

# Closed ranges where most contacts belong; to colour points differently
s.range = 'subj.range = c(-0.0001,0.004)'
r.range = 'ref.range = NULL' #c(-0.02,1)'

confMxMetric.v = c("MCC", "FDR", "TPR", "TNR", "PPV", "NPV")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#library(alphashape3d)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
source(paste0(wk.dir, "/lib/confusionMxMetric.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if( length(Metric$ref)!=1 ){ stop("Wrong Metric$ref.") }

ct.v.len <- length(ct.v)
m.len <- length(confMxMetric.v)
ct.v <- rep(ct.v, each=m.len)
confMxMetric.v <- rep(confMxMetric.v, times=ct.v.len)
x.len <- length(confMxMetric.v)

temp <- unlist(Metric)
id <- paste(paste(names(temp), temp, sep="_"), collapse="_")
run.name <- paste(src.id, id, sep="_"); rm(id, temp)
print(paste0(run.name, "..."), quote=FALSE)

out.name <- paste(gcb, chr, src.id, out.id, sep="_")
print(paste0("Output file name: ", out.name, "..."), quote=FALSE)

print(s.range, quote=FALSE)
print(r.range, quote=FALSE)

eval(parse(text=s.range))
eval(parse(text=r.range))

tbl <- data.matrix(expand.grid(1:length(Metric$subj), 1:length(Metric$ref)))
tbl.len <- nrow(tbl)
dimnames(tbl)[[2]] <- c("subj", "ref") 

Calc <- frSUBJ <- frREF <- list()
for(x in 1:x.len){
  
  metric <- confMxMetric.v[x]
  ct <- ct.v[x]
  csvnme <- paste(gcb, chr, ct, src.id, sep="_")
  
  temp <- list()
  for(i in 1:tbl.len){
    ind.v <- unname(tbl[i,])
    subj <- Metric$subj[ind.v[1]]
    ref <- Metric$ref[[ ind.v[2] ]]
    COMPIJMX <- read.csv(file=paste0(csv.dir, "/", csvnme, "_subj_", subj, "_ref_", 
                                     ref, ".csv"), header=TRUE, 
                         stringsAsFactors=FALSE)
    
    frSUBJ[[subj]] <- cbind.data.frame(id=subj, cutoff=COMPIJMX$c.offsubj, 
                                       fr=COMPIJMX$SP/COMPIJMX$final.nonNA.ij)
    frREF[[ct]] <- cbind.data.frame(id=ct, cutoff=COMPIJMX$c.offref,
                                    fr=COMPIJMX$RP/COMPIJMX$final.nonNA.ij)

    if(x==1){
      Calc[["cutoff"]] <- cbind(subj=COMPIJMX$c.offsubj, ref=COMPIJMX$c.offref)
    }
    id <- paste0(subj, "VS", ref)
    temp[[id]] <- confMxMetric(CONFMX=COMPIJMX[,c("TP","FP", "TN", "FN", "SP", "SN", 
                                                  "RP", "RN"),], metric=metric)
    print(paste0(id, " ", metric, " done!"), quote=FALSE)
    rm(COMPIJMX, ind.v, subj, ref, id)
  }
  
  Calc[[paste0(ct, "_", metric)]] <- do.call("cbind", temp)
  rm(temp, metric)
  
}

# Filter based on subj and ref cut-off ranges
if( is.null(subj.range) ){
  subj.range <- c( min(Calc$cutoff[,"subj"]), max(Calc$cutoff[,"subj"]) )
}
if( is.null(ref.range) ){
  ref.range <- c( min(Calc$cutoff[,"ref"]), max(Calc$cutoff[,"ref"]) )
}
out.range <- (Calc$cutoff[,"subj"]<subj.range[1] | Calc$cutoff[,"subj"]>subj.range[2]) |
  (Calc$cutoff[,"ref"]<ref.range[1] | Calc$cutoff[,"ref"]>ref.range[2])
num.incl <- sum(out.range==FALSE)
print(paste0(num.incl, "/", length(out.range), " rows after cut-off filtering..."), quote=FALSE)
Calc <- lapply(X=Calc, FUN=function(x){ x[out.range==FALSE, ] })

frSUBJ <- lapply(X=frSUBJ, FUN=function(x){ x[out.range==FALSE, ] })
frREF <- lapply(X=frREF, FUN=function(x){ x[out.range==FALSE, ] })

# Initialise output table
VOL.MX <- matrix(data=NA, nrow=length(confMxMetric.v), ncol=tbl.len)
dimnames(VOL.MX) <- list( names(Calc)[-1], dimnames(Calc[[2]])[[2]] )

# Per metric, choose rows of finite values
for( n in names(Calc)[-1] ){
  
  mx.TF <- apply(X=Calc[[n]], MARGIN=c(1,2), FUN=is.finite)
  incl.TF <- rowSums(mx.TF)==tbl.len
  print(paste0(n, ": ", sum(incl.TF), "/", num.incl, 
               " rows after finite values filtering..."), quote=FALSE)
  Calc[[n]] <- Calc[[n]][incl.TF,]
  
  VOL.MX[n,] <- colSums(Calc[[n]])
  
}

VOL.MX <- cbind.data.frame(ct=ct.v, metric=confMxMetric.v, VOL.MX, stringsAsFactors=FALSE)
max.ind <- apply(X=VOL.MX[,-(1:2)], MARGIN=1, FUN=function(rw){ which(max(rw)==rw) })
VOL.MX$max <- colnames(VOL.MX)[-(1:2)][max.ind]
write.csv(VOL.MX, file=paste0(out.dir, "/", out.name, ".csv"), row.names=TRUE)

# Plot fraction of contacts
frSUBJ <- do.call("rbind.data.frame", frSUBJ)
frREF <- do.call("rbind.data.frame", frREF)
coul0 <- rev(brewer.pal(n=11,name="Spectral"))
p.lst <- list()
for(x in c("SUBJ", "REF")){
  eval(parse(text=paste0(
    "p.lst[[x]] <- ggplot(data=fr", x, ", aes(x=cutoff, y=fr))"
  )))
  eval(parse(text=paste0(
    "coul <- colorRampPalette(coul0)(length(levels(fr", x, "$id)))"
  )))
  p.lst[[x]] <- p.lst[[x]] +
    geom_point(aes(col=id), alpha=0.5, size=2.5, shape=1) +
    scale_colour_manual(values=coul) +
    guides(color=guide_legend(ncol=1)) +
    labs(title=paste0(out.name, "_", x), 
         x="Cut-off", y="Fraction of contacts") +
    bgr2 + 
    theme(legend.text=element_text(size=15))
}
p.arr <- ggarrange(plotlist=p.lst, nrow=1, ncol=2)
ggexport(p.arr, height=10, width=20, filename=paste0(out.dir, "/", out.name, "_frCombined.pdf" ))

# rm(list=ls()); gc()

#x <- data.matrix(x[is.finite(x$densval),])
#dimnames(x)[[2]] <- c("", "", "Z")
#ashape3d.obj <- ashape3d(x=x, alpha=10, pert=TRUE)
#plot(ashape3d.obj, byComponents=FALSE, indexAlpha = 1)
#vol <- volume_ashape3d(ashape3d.obj, byComponents = FALSE, indexAlpha = 1)
