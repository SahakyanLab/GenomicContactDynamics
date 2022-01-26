################################################################################
# Make scatter plot of central value of gene lengths vs. Cp.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/16_GeneVsPersist"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lencp.dir = paste0(wk.dir, "/out_geneLengthPlot")
out.dir = paste0(wk.dir, "/out_geneLengthPlotLine")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
refseq = "ALL" 
out.id = "length" # "transcript"| "repeat"
#lentyp.v = c("len_full", "len_exon", "len_intron",
#             "mean_exon", "mean_intron", "div_intronBYexon",
#            "num_exon", "num_intron")
#plot.id = "transcript_lengths"

lentyp.v = c("len_repeat_full", "len_repeat_exon", "len_repeat_intron",
             "fr_repeat_full", "fr_repeat_full", "fr_repeat_full",
             "len_full", "len_exon", "len_intron")
plot.id = "transcript_repeat_values"

repfree = F

centr = "mean" # "median"

col.v = c(full="#FDC776", exon="#e25d6c", intron="#3288BD", intronBYexon="#bb88dd")
# len=15
shp.v = c(`fr_repeat`=13, `len_repeat`=11, len=11, num=0, mean=8, div=17)
facet = T
ylim.v = c(-0.1,2)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
id <- paste0("hg19anno_", refseq, "_", gcb)
out.id <- paste0(id, "_", plot.id, "_", centr, "_repfree", repfree, "_facet", facet)

DF <- list()
for(lentyp in lentyp.v){
  
  load(file=paste0(lencp.dir, "/", id, "_", lentyp, ".RData"))
  cp.v <- sort(unique(LENCP.DF[,"cp"]), decreasing=F)
  
  LENCP.DF <- as.data.frame(LENCP.DF)
  
  LENCP.DF$cp <- factor(x=as.character(LENCP.DF$cp), 
                        levels=as.character(cp.v))
  
  var.v <- colnames(LENCP.DF)[colnames(LENCP.DF)!="cp"]
  
  if(repfree==F){
    var.v <- var.v[!grepl(x=var.v, pattern="_repfree", fixed=T)]
  }
  
  for(vr in var.v){
    
    data.table::setnames(x=as.data.frame(LENCP.DF), old=vr, new="L", 
                         skip_absent=F)
    eval(parse(text=paste0(
      paste0( 'tmp <- aggregate(formula=L~cp, data=LENCP.DF, FUN=', centr, ')' )
    )))
    
    tmp$L.center <- tmp$L 
    tmp$L <- log2(tmp$L/tmp[as.character(tmp$cp)=="1","L"])
    DF[[vr]] <- cbind(L.name=rep(vr), tmp)
    
    LENCP.DF$L <- NULL
    
    rm(vr, tmp)
    
  }
  
  print(paste0(lentyp, " done!"), quote=FALSE)
  
  rm(LENCP.DF, cp.v, var.v)
  
} # length.v for loop end

DF <- do.call("rbind", DF)
rownames(DF) <- NULL

DF$L.name <- factor(x=DF$L.name, levels=unique(DF$L.name))
grp <- as.character(DF$L.name)

DF$repcontgrp <- "orig"
DF$repcontgrp[grepl(x=DF$L.name, pattern="_repfree", fixed=T)] <- "repfree"

# Groups are div, len, len_repeat, mean, num, fr_repeat
calcgrp <- strsplit(x=grp, split="_full|_exon|_intron", fixed=F)
calcgrp <- unlist(lapply(X=calcgrp, FUN=function(x) x[1]))
DF$calcgrp <- calcgrp
DF$calcgrp <- factor(x=DF$calcgrp, 
                     levels=intersect(names(shp.v), unique(DF$calcgrp)))

# Groups are exon, full, intron, intronBYexon
partgrp <- mapply(FUN=gsub, x=grp, pattern=paste0(calcgrp, "_|_repfree"), 
                  replacement="", fixed=F, SIMPLIFY=T)
DF$partgrp <- partgrp
DF$partgrp <- factor(x=DF$partgrp, 
                     levels=intersect(names(col.v), unique(DF$partgrp)))

rm(grp, calcgrp, partgrp)

# Plot
p <- ggplot(data=DF, aes(x=cp, y=L, group=L.name)) +
  geom_line(colour="gray91", size=1) +
  geom_point(aes(shape=calcgrp, colour=partgrp), size=4) +
  geom_hline(yintercept=0, linetype="dashed", colour="gray70", size=0.5) + 
  scale_y_continuous(limits=ylim.v) + 
  scale_colour_manual(values=col.v[levels(DF$partgrp)]) +
  scale_shape_manual(values=shp.v[levels(DF$calcgrp)]) + 
  labs(title=out.id, x="cp", colour="",
       y=paste0("log2 or FC ", centr, " wrt cp=1")) +
  bgr2 +
  theme(legend.text=element_text(size=15, face="bold"),
        legend.title=element_text(size=20, face="bold"),
        legend.key = element_rect(colour="transparent",
                                  fill="transparent"),
        strip.background=element_blank())
   
if(facet){
  p <- p + facet_grid(.~repcontgrp)
}

ggsave(filename=paste0(out.dir, "/", out.id, "_lineplot.pdf"),
       width=15, height=10, units="in", plot=p)

# rm(list=ls()); gc()