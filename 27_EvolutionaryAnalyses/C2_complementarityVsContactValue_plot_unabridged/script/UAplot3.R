################################################################################
# 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
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
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/27_EvolutionaryAnalyses")
valmx.dir = paste0(wk.dir, "/out_complementarityVsContactValue_unabridged") 
out.dir = paste0(wk.dir, "/out_complementarityVsContactValue_plot_unabridged") 
state.group.file = paste0(wk.dir, "/Phylo-HMRF_state_group.csv")
group.col.file = paste0(wk.dir, "/Phylo-HMRF_state_group_col.csv")
### OTHER SETTINGS #############################################################
value.id = "min2Mb_genome_state_Phylo-HMRF_contact50K_norm_KR_CPhg38ToHg19_LOwidth.min.bp30000_kmer_min2Mb_consensusCp_CFcontact50K_NONE_KR"
state.group = "Grp4"
consensus.Cp.val = "min.consCp" #"median.consCp"
chrs = paste0("chr", c(1:22, "X"))
bin.len = 50000
gap.range.bins.closed = NULL # j - i - 1, No gap filtering if NULL
out.id = paste0("hg38.108_bin", bin.len, "bp_", value.id, "_gaprangebins_", 
                gap.range.bins.closed[1], "_", gap.range.bins.closed[2],
                "_state", state.group, "_XEQUALSceilingOF", consensus.Cp.val)  
grp.order.bar = rev(c("C-high", "C-mid", "C-low", "WC", "NC-LS", "NC-hom_high", 
                      "NC-hom_low", "NC", "NONE"))
coul.spec <- c(`GSM3685694_panTro5_mapping`="cyan2", 
               `GSM3685695_panPan2_mapping`="blue", 
               `GSM3685696_gorGor4_mapping`="darkgreen", 
               `hg38.downsample`="hotpink") 
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(RColorBrewer)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/compareManyDist.R"))
### FUNCTION ###################################################################
countPerX <- function(a){
  a <- a[!is.na(a)]
  return( c(y=min(a), label=length(a)) )
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
df <- sapply(chrs, simplify=F, FUN=function(chr){
  
  load(paste0(valmx.dir, "/", chr, "_", value.id, ".RData"))
 
  if( !is.null(gap.range.bins.closed) ){
    print(paste0(chr, ": Gap Filtering"))
    gaps <- unname(VAL.MX$j - VAL.MX$i - 1)
    VAL.MX <- VAL.MX[gaps >= gap.range.bins.closed[1] & gaps <= gap.range.bins.closed[2],]
  }
  
  print(paste0(chr, ": Data obtained."))
  return(VAL.MX)
  
})
df <- do.call("rbind.data.frame", df)

meta <- lapply(df, FUN=function(col){
  sum(is.finite(col))
})
meta <- stack(meta)
colnames(meta) <- c("countFiniteRow", "column name")
write.table(meta, file=paste0(out.dir, "/", out.id, "_countPerColumnInDF_metadata.txt"),
            row.names=F, quote=F, sep="\t")

# Change "||" to "II" in column names for ggplot
colnames(df) <- gsub(pattern="||", replacement="II", x=colnames(df), fixed=T)

# Identify PHYLO-HMRF state groupings 
grp.df <- read.csv(state.group.file, header=T, stringsAsFactors=F)
is.nonNAstate <- !is.na(df$PHYLOHMRFstate)
df$grp[is.nonNAstate] <- grp.df[[state.group]][ as.numeric(as.character(df$PHYLOHMRFstate))[is.nonNAstate] ]

# Get Cps for plotting
df$Cp.plot <- ceiling(df[[consensus.Cp.val]])
df$Cp.plot <- factor(as.character(df$Cp.plot), 
                     levels=as.character(sort(unique(df$Cp.plot))))

## PLOTS

P.LST <- list()

# Group/state colour
col.df <- read.csv(group.col.file, header=T, stringsAsFactors=F)
col.df$Col <- paste0("#", col.df$Col)

# Violin plot: Cp vs. complementarity values in hg38 genome

compl.nmes <- colnames(df)
compl.nmes <- compl.nmes[ grepl(pattern="kmer.|align.", compl.nmes) ]
for(nme in compl.nmes){
  
  df.tmp <- na.omit(df[,c("Cp.plot", nme)])
  compareManyDist(xval=df.tmp[[nme]], grp=df.tmp$Cp.plot, 
                  alt="two.sided", out.dir=out.dir, out.name=paste0(out.id, "_", nme))

  P.LST[[nme]] <- ggplot(data=df.tmp, aes_string(x="Cp.plot", y=nme)) +
    geom_violin(scale="width", fill="#FDC776", col="#FDC776", trim=T, lwd=2) +
    stat_boxplot(geom="errorbar", width=0, lty="dashed") + 
    stat_summary(fun="median", pch=1, position=position_dodge(0.5)) +
    labs(x=NULL, title=paste0(out.id, 
                              "\n geom_violin(trim=T) so tails at range of data_hollow point at median_vertical line is default stat_boxplot errorbar that is same as to vertical line when making a default ggplot2 boxplot")) + 
    bgr2 + 
    theme(legend.position="none", plot.title=element_text(size=7),
          axis.text.x=element_text(size=7),
          aspect.ratio=0.25) 
  P.LST[[paste0(nme, "withCount")]] <- P.LST[[nme]] + 
    stat_summary(fun.data=countPerX, geom="text") 
  
  rm(df.tmp)
  
}

# Median + tails plot: Cp vs. contact frequency data

KR.nmes <- colnames(df)[grepl(colnames(df), pattern="_KR", fixed=T)] 
KR.df <- df[,c(KR.nmes, "Cp.plot")]
KR.df <- reshape2::melt(KR.df, id="Cp.plot")

coul.nme <- gsub(pattern="_KR", replacement="", x=levels(KR.df$variable))
P.LST[["KR_Cp"]] <- ggplot(data=KR.df, aes(x=Cp.plot, y=log10(value))) +
  stat_boxplot(aes(col=variable), geom="errorbar", width=0, lty="dashed", position=position_dodge(0.5)) +
  stat_summary(aes(col=variable), fun="median", pch=1, position=position_dodge(0.5)) +
  scale_x_discrete(limits=levels(KR.df$Cp.plot)) +
  scale_colour_manual(values=unname(coul.spec[coul.nme])) +
  labs(x=NULL, y="log10(KR)", 
       title=paste0(out.id, 
                    "\n note that KR/NONE cf values not cross-species comparable_hollow point at median_vertical line is default stat_boxplot errorbar that is same as to vertical line when making a default ggplot2 boxplot")) +
  bgr2 +
  theme(legend.position="bottom", legend.text=element_text(size=2), 
        legend.title=element_text(size=2), plot.title=element_text(size=7),
        axis.text.x=element_text(size=7), aspect.ratio=0.25) 

P.LST[[paste0("KR_Cp_withCount")]] <- P.LST[["KR_Cp"]] + 
  stat_summary(fun.data=countPerX, geom="text")

for(nme in KR.nmes){
  
  df.tmp <- na.omit(df[,c("Cp.plot", nme)])
  compareManyDist(xval=df.tmp[[nme]], grp=df.tmp$Cp.plot, 
                  alt="two.sided", out.dir=out.dir, out.name=paste0(out.id, "_", nme))
  rm(df.tmp)
    
}

rm(KR.df, coul.nme)

# Bar-plot: Phylo-HMRF states vs. Cp

bar.df <- df[,c("grp", "Cp.plot")]
bar.df$grp[is.na(bar.df$grp)] <- "NONE"

bar.df <- na.omit(bar.df) 
bar.df$grp <- factor(bar.df$grp, 
                     levels=intersect(grp.order.bar, unique(bar.df$grp)))

barcol.df <- col.df[col.df$Grp %in% levels(bar.df$grp),]
barcoul.grp <- setNames(barcol.df$Col, nm=barcol.df$Grp)

P.LST[["grp_Cp_bar"]] <- ggplot(data=bar.df, aes(x=Cp.plot)) +
  geom_bar(aes(fill=grp), position="fill") + 
  scale_fill_manual(values=unname(barcoul.grp[levels(bar.df$grp)])) + 
  labs(x=NULL, title=out.id) + 
  guides(fill=guide_legend(ncol=length(levels(bar.df$grp)))) + 
  bgr2 + 
  theme(legend.position="bottom", legend.text=element_text(size=2), 
        legend.title=element_text(size=2), plot.title=element_text(size=7),
        axis.text.x=element_text(size=7), aspect.ratio=0.25) 

P.LST[["grp_Cp_bar_withCount"]] <- ggplot(data=bar.df, aes(x=Cp.plot)) +
  geom_bar() +
  labs(title="grp_Cp_bar_withCount") +
  geom_text(stat="count", label=unname(table(bar.df$Cp.plot)), size=3, vjust=-0.2) +
  bgr2 + 
  theme(legend.position="bottom", legend.text=element_text(size=2), 
        legend.title=element_text(size=2), plot.title=element_text(size=7),
        axis.text.x=element_text(size=7), aspect.ratio=0.25) 

rm(bar.df, barcol.df, barcoul.grp)

# Bar-plot: across-species Cp vs. Cp

NONE.nmes <- colnames(df)[grepl(colnames(df), pattern="_NONE", fixed=T)] 

df$Cp.spec.KR <- unname( rowSums(!is.na(df[,KR.nmes])) )
df$Cp.spec.NONE <- unname( rowSums(!is.na(df[,NONE.nmes])) )

# Comparison plot of across species Cp based on KR vs. NONE
diff.val <- diff(df$Cp.spec.NONE-df$Cp.spec.KR)
diff.meta <- stack(table(diff.val))
colnames(diff.meta) <- c("count", "dfCp.spec.NONEMINUSdfCp.spec.KR")
write.table(diff.meta, file=paste0(out.dir, "/", out.id, "_dfCp.spec.NONEMINUSdfCp.spec.KR__metadata.txt"),
            row.names=F, quote=F, sep="\t")
pdf(paste0(out.dir, "/", out.id, "_densityOfdfCp.spec.NONEMINUSdfCp.spec.KR.pdf"), 
    height=10, width=10)
plot(density(diff.val), main="diff(df$Cp.spec.NONE-df$Cp.spec.KR)")
dev.off()

cpspec.df <- df[,c("Cp.plot", "Cp.spec.KR", "Cp.spec.NONE")]
cpspec.df$Cp.spec.KR <- factor(as.character(cpspec.df$Cp.spec.KR), 
                               levels=as.character(sort(unique(cpspec.df$Cp.spec.KR))))
cpspec.df$Cp.spec.NONE <- factor(as.character(cpspec.df$Cp.spec.NONE), 
                                 levels=as.character(sort(unique(cpspec.df$Cp.spec.NONE))))
coul <- colorRampPalette( rev(brewer.pal(11, "Spectral")) )(4)
coul <- c("black", coul)  
names(coul) <- levels(cpspec.df$Cp.spec.KR)

for( cf.nme in c("KR", "NONE") ){
  
  col.nme <- paste0("Cp.spec.", cf.nme)
  df.tmp <- cpspec.df[,c("Cp.plot", col.nme)]
  df.tmp <- na.omit(df.tmp)
  P.LST[[paste0(col.nme, "_Cp_bar")]] <- ggplot(data=df.tmp, aes(x=Cp.plot)) +
    geom_bar(aes_string(fill=col.nme), position="fill") + 
    scale_fill_manual(values=unname( coul[levels(df.tmp[[col.nme]])] )) + 
    labs(x=NULL, title=out.id, fill=cf.nme) + 
    guides(fill=guide_legend( ncol=length(levels(df.tmp[[col.nme]])) )) + 
    bgr2 + 
    theme(legend.position="bottom", legend.text=element_text(size=2), 
          legend.title=element_text(size=2), plot.title=element_text(size=7),
          axis.text.x=element_text(size=7), aspect.ratio=0.25) 
  
  P.LST[[paste0(col.nme, "_Cp_bar_withCount")]] <- ggplot(data=df.tmp, aes(x=Cp.plot)) +
    geom_bar() +
    labs(title=paste0(col.nme, "_Cp_bar_withCount")) +
    geom_text(stat="count", label=unname(table(df.tmp$Cp.plot)), size=3, vjust=-0.2) +
    bgr2 + 
    theme(legend.position="bottom", legend.text=element_text(size=2), 
          legend.title=element_text(size=2), plot.title=element_text(size=7),
          axis.text.x=element_text(size=7), aspect.ratio=0.25) 
  
  rm(df.tmp)
  
}

# Save plots

for( p.nme in names(P.LST) ){
  
  print(paste0("Plotting ", p.nme, "..."), quote=F)
  ggsave(paste0(out.dir, "/", out.id, "_", p.nme, ".pdf"), plot=P.LST[[p.nme]], 
         width=15, height=5, units="in")

}

# rm(list=ls()); gc()
