################################################################################
# 
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
src.dir = out.dir = paste0(wk.dir, "/z_ignore_git/out_combine")
#elm.file = paste0(src.dir, "/out_filterElmTRUE_repName.txt")
ijcount.file = paste0(wk.dir, "/min2Mb_ij_2Mb_Cp1To21_ijcount.txt")
### OTHER SETTINGS #############################################################
ijcount.Cps = 1:21
est.src.id = "metric_estimatemedian_Cpdyn1T3_per19To21"
withLabel = T
out.affix = paste0("withLabel", withLabel)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(reshape2)
library(cowplot)
library(ggrepel)
#library(RColorBrewer)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
#theme.nolabels <- list(axis.title.y=element_blank(), axis.text.x=element_blank(),
#                       axis.text.y=element_blank())

#elements <- readLines(elm.file)
load(paste0(src.dir, "/repeat_group.RData"))
elements <- df$repName[ !df$plot.group %in% "not.transposon" ]
#elements <- unique(df$repName)
rm(df)
elm.order <- elements # Order of elements in estimate plots

# Count of contacts per Cp
ij.count <- setNames(as.numeric(readLines(ijcount.file)), nm=ijcount.Cps)

##### Proportion of contacts with at least 2 sites across Cp, plotted as fold change
##### relative to Cp=1 value

load(paste0(src.dir, "/minrep_contact_coverage.RData"))
mx <- mx[elements,]
mx <- mx[order(mx[,"21"], decreasing=T),]

# Fraction of contacts per Cp (use transpose if mx is not big)
if( identical(dimnames(mx)[[2]], names(ij.count)) ){
  
  fr.mx <- t( t(mx) / ij.count )
  
  #fr.mx <- mx
  #for(cp in as.character(ijcount.Cps)){
  #  fr.mx[,cp] <- mx[,cp] / ij.count[[cp]]
  #}
  
} else {
  stop("ijcov: Wrong element order.")
}

# Calculate fold change relative to Cp=1 value

fc.mx <- log2(fr.mx / fr.mx[,"1"])
fc.df <- reshape2::melt(fc.mx)
setnames(fc.df, old=c("Var1", "Var2"), new=c("repName", "Cp"))

# Add proportion at Cp=21 as label

fr.df <- reshape2::melt(fr.mx)
setnames(fr.df, old=c("Var1", "Var2"), new=c("repName", "Cp"))

if( identical(fc.df[,1:2], fr.df[,1:2]) ){
  
  fc.df$label <- NA
  is.cp21 <-fc.df$Cp == 21
  fc.df$label[is.cp21] <- format(fr.df$value[is.cp21], digits=2)
  
} else {
  stop("ijcov: fc.df and fr.df factors not identical.")
}

# Plot

fc.df$repName <- factor(fc.df$repName, levels=rownames(mx))
fc.df$Cp <- factor(as.character(fc.df$Cp), levels=ijcount.Cps)

p <- ggplot(data=fc.df, aes(x=Cp, y=value, group=repName)) +
  geom_hline(yintercept=0, linetype="dashed") + 
  geom_line(aes(col=repName)) + 
  geom_point(aes(col=repName)) +
  #geom_text_repel(aes(label=label), size=3.5, box.padding=0.5) + 
  #scale_colour_brewer(palette="Paired") + 
  labs(y="fc", 
       title="minrep_contact_coverage_labelFractionAtCp=21_SubfamInLegendInDecreasingCp=21Fraction") + 
  bgr2 +
  theme(plot.title=element_text(size=10))

ggsave(filename=paste0(out.dir, "/minrep_contact_coverage_fc_", out.affix, ".pdf"),
       plot=p, height=10, width=10)

##### Scatter comparing ijcov Cp19To21 with Cp1To3

load(paste0(src.dir, "/minrep_contact_coverage.RData"))
mx <- mx[elements,]
mx <- mx[order(mx[,"21"], decreasing=T),]

mx <- cbind( dyn=rowSums(mx[,as.character(1:3)]), 
             per=rowSums(mx[,as.character(19:21)]) )
mx[,"dyn"] <- mx[,"dyn"] / sum(ij.count[as.character(1:3)])
mx[,"per"] <- mx[,"per"] / sum(ij.count[as.character(19:21)])

DF <- as.data.frame(mx)
rm(mx)
DF$repName <- rownames(DF)
#DF$repName <- factor(DF$repName, levels=rownames(mx))
DF$plot.group <- rep("ijcov", times=nrow(DF))
rownames(DF) <- NULL

##### Scatter comparing minrep.nonzero.fr Cp19To21 with Cp1To3
load(paste0(src.dir, "/metric_nonzero_Cpdyn1T3_per19To21.RData"))
DF1 <- as.data.frame(mx)
rm(mx)
DF1 <- DF1[elm.order,]
colnames(DF1) <- substr(colnames(DF1), start=1, stop=3)

DF1$dyn <- DF1$dyn / sum(ij.count[as.character(1:3)])
DF1$per <- DF1$per / sum(ij.count[as.character(19:21)])

DF1$repName <- elm.order
DF1$plot.group <- "aminrepnon0.fr"
rownames(DF1) <- NULL

##### Scatter comparing metric estimate Cp19To21 with Cp1To3

metrics <- c("min", "skew", "sum")
metric.len <- length(metrics)
DF2 <- list()
for(m in 1:metric.len){
  
  metric.id <- paste0(metrics[[m]], "rep")
  load(paste0(src.dir, "/", metric.id, "_", est.src.id, ".RData"))
  mx <- mx[elm.order,]
  dimnames(mx)[[2]] <- substr(colnames(mx), start=1, stop=3)
  DF2[[m]] <- as.data.frame(mx)
  DF2[[m]]$repName <- dimnames(mx)[[1]]
  DF2[[m]]$plot.group <- rep(paste0(metric.id, ".est"), 
                             times=nrow(DF2[[m]]))

  rm(mx)
  
}

DF2 <- do.call("rbind", DF2)
rownames(DF2) <- NULL

DF <- rbind(DF, DF1, DF2)
# DF$subfam.grp <- NA
# DF$subfam.grp[grepl("MIR", DF$repName)] <- "MIR" 
# DF$subfam.grp[grepl("Alu", DF$repName)] <- "Alu" 
# DF$subfam.grp[grepl("L2", DF$repName)] <- "L2" 

##### Plot

DF$repName <- factor(DF$repName, levels=elm.order)

# elm.cols <- c(brewer.pal(n=9,"Purples")[c(3,6,8)],
#               brewer.pal(n=9,"Greens")[c(8,3,6)],
#               brewer.pal(n=9,"Reds")[c(6,2,8,4)]
# )

# elm.cols <- c(brewer.pal(n=9,"Purples")[c(8,8,8)],
#                brewer.pal(n=9,"Greens")[c(8,8,8)],
#                brewer.pal(n=9,"Reds")[c(8,8,8,8)]
# )

## Plot contact coverage

df.plot <- DF[DF$plot.group %in% c("ijcov", "aminrepnon0.fr"),]
df.plot$label <- NA
is.label <- !is.na(df.plot$per) & df.plot$per > 0.1
df.plot$label[is.label] <- as.character(df.plot$repName)[is.label]
p <- ggplot(data=df.plot, aes(x=dyn, y=per)) +
  geom_hline(yintercept=0.5, linetype="dashed") +
  geom_vline(xintercept=0.5, linetype="dashed") +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="black") +
  #geom_hex() +
  #stat_binhex(aes(label=..count..), geom="text", colour="white") +
  geom_point(shape=20, size=10, col=adjustcolor("black", 0.1)) +
  geom_point(shape=20, size=2, col=adjustcolor("black", 0.5)) +
  scale_x_continuous(limits=c(0,1.0)) +
  scale_y_continuous(limits=c(0,1.0)) +
  #scale_fill_gradient(low="#3C5488FF", high="#DC0000FF") + 
  bgr2 +
  facet_grid(.~plot.group)

if(withLabel){
  p <- p +
    geom_text_repel(aes(label=label), box.padding=1, col="darkgreen",
                    size=1, segment.size=1, min.segment.length=0) 
} else {
  p <- p + 
    theme(axis.title.y=element_blank(), axis.text.x=element_blank(), 
          axis.text.y=element_blank())
}

ggsave(filename=paste0(out.dir, "/scatter_dynCp1To3VsPerCp19To21_", out.affix, ".png"),
       plot=p, height=10*300, width=20*300, units="px")

rm(df.plot)

#

# elm.cols <- c(brewer.pal(n=9,"Purples")[c(8,8,8)],
#               brewer.pal(n=9,"Greens")[c(8,8,8)],
#               brewer.pal(n=9,"Reds")[c(8,8,8,8)]
# )

## Plot metric estimates

DF$repName <- as.character(DF$repName)
DF$repName <- factor(DF$repName, levels=sort(unique(DF$repName)))

metrics <- c("minrep.est", "skewrep.est", "sumrep.est")
axis.lims <- list(c(-0.25,3), c(0,1), c(-1,6))
P.LST <- sapply(X=1:metric.len, simplify=F, FUN=function(m){
  
  metric <- metrics[[m]]
  
  p <- ggplot(data=DF[DF$plot.group %in% metric,], aes(x=dyn, y=per)) +
    geom_abline(slope=1, intercept=0, linetype="dashed") +
    #geom_point(shape=20, size=10) +
    #geom_hex() +
    #stat_binhex(aes(label=..count..), geom="text", colour="white") +
    geom_point(shape=20, size=10, col=adjustcolor("black", 0.1)) +
    geom_point(shape=20, size=2, col=adjustcolor("black", 0.5)) +
    scale_x_continuous(limits=axis.lims[[m]]) +
    scale_y_continuous(limits=axis.lims[[m]]) +
    scale_shape_manual(values=c(1,3,6)) + 
    #scale_fill_gradient(low="#3C5488FF", high="#DC0000FF") + 
    labs(y=paste0(metric, "_perVsdyn")) +
    bgr2 + 
    theme(legend.position="none", plot.title=element_blank(), 
          axis.title.x=element_blank(), axis.title.y=element_text(size=10))
  
  if(withLabel){
    p <- p + 
      geom_text_repel(aes(label=repName), box.padding=1, col="darkgreen",
                      size=3, segment.size=1, min.segment.length=1) 
  } else {
    p <- p + 
      theme(axis.title.y=element_blank(), axis.text.x=element_blank(), 
            axis.text.y=element_blank())
  }
    
  return(p)
  
})

row.num <- 2
col.num <- 2
p.arr <- plot_grid(plotlist=P.LST, nrow=row.num, ncol=col.num, align="hv", byrow=T)
save_plot(filename=paste0(out.dir, "/scatter_metric_estimates_dynCp1To3VsPerCp19To21_", out.affix, ".png"), 
          plot=p.arr, base_width=col.num * 10, base_height=row.num * 10)


# > COPYNUM$repName[is.elm,]
# repName copyNumber   fraction
# 1111    MIRb     223577 0.04273067
# 1109     MIR     174175 0.03328882
# 645      L2a     169322 0.03236130
# 319    AluJb     142591 0.02725240
# 334    AluSx     141949 0.02712970
# 647      L2c     139481 0.02665800
# 340     AluY     118506 0.02264920
# 1112    MIRc     102688 0.01962602
# 646      L2b      97239 0.01858459
# 1110    MIR3      90185 0.01723641

# rm(list=ls()); gc()