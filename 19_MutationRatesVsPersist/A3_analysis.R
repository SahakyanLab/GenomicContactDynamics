################################################################################
# Analyse trends per signature
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
plot.dir = paste0(wk.dir, "/out_plotsVsCp")
out.dir = paste0(wk.dir, "/out_analysis")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
data.id = "donor_centric_PCAWG" # "CosmicNCV", "donor_centric_PCAWG"
src.id = "Hg19" # "Hg19" | "hg38ToHg19"
mut.v = c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
calc = "Nmsitenorm" #c("Tmut", "Nmsite", "TmutDIVNmsite", "Nmsitenorm", "numWTSEQ")
aggregatefunx = "median"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(ggpubr)
library(ggsci)
library(ComplexHeatmap)
library(gridGraphics)
library(grid)
library(gridExtra)
source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
makeScatter <- function(df, calc, yint=NULL, plot.id){
  
  # Plot
  df$pval <- c(`0`="ns", `1`="sig")[ as.character(as.numeric(df$pval<0.05)) ]
  df$pval <- factor(df$pval, levels=c("ns", "sig"))
  alpha.v <- c(0.3, 1)
  names(alpha.v) <- levels(df$pval)
  
  if( !is.character(df$ind) ){
    stop("Scatterplot making: ind not a character as expected.")
  }
  df$ind <- factor(as.character(df$ind), 
                   levels=c(as.character(sort(as.numeric(unique(df$ind)))))
  )
  
  df$alt <- factor(as.character(df$alt), 
                   levels=c("greater", "less", "two.sided"))
  
  df$loc.id <- factor(as.character(df$loc.id),
                      levels=c("exon", "intron", "intergenic", "intron_intergenic", 
                               "exon_intron_intergenic_intergenic.excl"))
  col.v <- ggsci::pal_npg("nrc")( length(levels(df$loc.id)) )
  names(col.v) <- levels(df$loc.id)
  
  #df$sigEpLim.id <- factor(as.character(df$sigEpLim.id), 
  #                         levels=c("nosampfilter", "-1_0", "10_100", "30_100", "50_100")
  #)
  shape.v <- 15:21
  #names(shape.v) <- levels(df$sigEpLim.id)
  
  p <- ggplot(data=df, aes(x=ind, y=values)) +
    geom_point(aes(colour=loc.id, alpha=pval, shape=sigEpLim.id), size=3) +
    scale_colour_manual(values=col.v) + 
    scale_alpha_manual(values=alpha.v) + 
    scale_shape_manual(values=shape.v) + 
    labs(x="Cp", y=calc, title=plot.id) + 
    bgr2 + 
    theme(plot.title=element_text(size=15), legend.text=element_text(size=10),
          legend.title=element_text(size=10))
  
  if( !is.null(yint) ){
    p <- p + geom_hline(yintercept=yint, linetype="dashed", colour="gray70", size=0.5) 
  }
  
  return(p)
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(gcb, "_", data.id, "_", src.id, "_", calc, "_", aggregatefunx)

mut.id.v <- sapply(X=mut.v, simplify=T, FUN=function(mut){
  gsub(x=mut, pattern=">", replacement="To", fixed=T)
})
mut.id.v.len <- length(mut.id.v)

DF <- list()
for( mut.id in unname(mut.id.v) ){
  
  #for(calc in  calc.v){
    
    pdata.id <- paste0(gcb, "_", data.id, "_", src.id, "_", mut.id, "_", calc, "_",
                       aggregatefunx)
    # P.DF
    load(file=paste0(plot.dir, "/", pdata.id, "_scatplot.RData"))
    # PFC.DF
    load(file=paste0(plot.dir, "/", pdata.id, "_FC_scatplot.RData"))
    
    DF[[paste0(mut.id, calc)]] <- cbind.data.frame(P.DF, fc=PFC.DF$values,
                                                   stringsAsFactors=F)
    
    rm(P.DF, PFC.DF); gc()
    
  #} # calc.v for loop end
  
} # mut.id.v for loop end

DF <- do.call("rbind.data.frame", c(DF, stringsAsFactors=F))
rownames(DF) <- NULL

# Remove duplicated nosampfilter entries
tmp <- DF[DF$sigEpLim.id=="nosampfilter",]
tmp <- tmp[!duplicated(tmp[,-2]),]
tmp$SIG.id <- "nosampfilter"

DF <- rbind.data.frame(tmp, DF[DF$sigEpLim.id!="nosampfilter",])
Cp.v <- as.character(sort(unique(DF$ind)))

# Plot trend for nosampfilter (all mutations)
tmp$ind <- as.character(tmp$ind)
tmp$sigEpLim.id <- tmp$mut.id
tmp$values <- tmp$fc

p <- makeScatter(df=tmp, calc=calc, yint=NULL, plot.id=out.id)
ggsave(filename=paste0(out.dir, "/", out.id, "_nosampfilterTrend.pdf"), 
       height=10, width=10, units="in")
rm(tmp, p); gc()

# Add trend column
DF$trend <- DF$alt
DF$trend[DF$pval>0.05]  <- "ns"
DF$trend[DF$trend=="greater"] <- "g"
DF$trend[DF$trend=="less"] <- "l"
#> table(DF$trend, useNA="always")
#g     l    ns  <NA> 
#  41086  3519 62495     0 

if( any(DF$trend=="two.sided" & is.na(DF$trend)) ){
  stop("Checkpoint 1.")
}

#-------------------percbinmut

loc.id.v <- unique(DF$loc.id)
loc.id.v.len <- length(loc.id.v)

pdf(file=paste0(out.dir, "/", out.id, "_percbinmutbp.pdf"), height=10,
    width=20)
par( mfrow=c(2,3) )
for(loc.id in loc.id.v){
  
  boxplot(percbinmut~sigEpLim.id, outline=FALSE, data=DF[DF$loc.id==loc.id,],
          xlab="sigEpLim.id", ylab="percbinmut", boxwex=0.6, cex.axis=1.2, 
          col="#FDC776", cex.main=1, main=paste0(out.id, "_", loc.id))
  
}
dev.off()

p <- ggplot(data=DF, aes(x=sigEpLim.id, y=percbinmut, fill=loc.id)) +
  geom_violin() +
  geom_boxplot(width=0.1, color="grey", alpha=0.2, outlier.colour="darkred") +
  scale_y_continuous(breaks=seq(0,100,20)) +
  scale_fill_npg() + 
  bgr2 +
  theme(axis.text.x=element_text(size=10)) +
  facet_wrap(.~loc.id) 
  
ggsave(filename=paste0(out.dir, "/", out.id, "_percbinmutvp.pdf"), height=10,
       width=20, units="in")

#-------------------General trend per mut.id, loc.id and sigEpLim.id

fill.v <- ggsci::pal_npg("nrc")(3)
names(fill.v) <- c("ns", "l", "g")
for(loc.id in loc.id.v){
  
  p.lst <-list()
  for(mut.id in unname(mut.id.v)){
    
    df <- DF[DF$loc.id==loc.id & DF$mut.id==mut.id,]
    df$ind <- factor(x=as.character(df$ind), levels=Cp.v)
    df$trend <- factor(x=as.character(df$trend), levels=c("ns", "l", "g"))
    
    p.lst[[mut.id]] <- ggplot(data=df, aes(fill=trend,  x=ind)) + 
      geom_bar(position="dodge") +
      scale_fill_manual(values=fill.v) + 
      labs(x="Cp", title=paste0(out.id,  "_", mut.id, "_", loc.id)) + 
      bgr2 +
      theme(panel.grid.major.x=element_line(size=0.5, colour="gray70",  
                                            linetype="dashed"),
            axis.text.x=element_text(size=5)) +
      facet_grid(cols=vars(sigEpLim.id))
  
    rm(df); gc()
    
  }
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=mut.id.v.len, ncol=1, legend=NULL)
  ggexport(p.arr, height=5*mut.id.v.len, width=15, 
           filename=paste0(out.dir, "/", out.id, "_", loc.id, "_trend.pdf"))
  
}

for( mut.id in unname(mut.id.v) ){
  
  p.lst <-list()
  for(loc.id in loc.id.v){
    
    df <- DF[DF$loc.id==loc.id & DF$mut.id==mut.id,]
    df$ind <- factor(x=as.character(df$ind), levels=Cp.v)
    df$trend <- factor(x=as.character(df$trend), levels=c("ns", "l", "g"))
    
    p.lst[[loc.id]] <- ggplot(data=df, aes(fill=trend,  x=ind)) + 
      geom_bar(position="dodge") +
      scale_fill_manual(values=fill.v) + 
      labs(x="Cp", title=paste0(out.id,  "_", mut.id, "_", loc.id)) + 
      bgr2 +
      theme(panel.grid.major.x=element_line(size=0.5, colour="gray70",  
                                            linetype="dashed"),
            axis.text.x=element_text(size=5)) +
      facet_grid(cols=vars(sigEpLim.id))
    
    rm(df); gc()
    
  }
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=loc.id.v.len, ncol=1, legend=NULL)
  ggexport(p.arr, height=5*loc.id.v.len, width=15, 
           filename=paste0(out.dir, "/", out.id, "_", mut.id, "_trend.pdf"))
  
}

#-------------------Trends per signature, heatmaps
DF <- DF[DF$SIG.id!="nosampfilter",]
sigEpLim.id.v <-unique(DF$sigEpLim.id)

if( any(!DF$trend%in%c("ns", "g", "l")) ){
  stop("Checkpoint 2.")
}
DF$trendnum <- 0
DF$trendnum[DF$trend=="g"] <- 1
DF$trendnum[DF$trend=="l"] <- -1

grab_grob <- function(){
  gridGraphics::grid.echo()
  grid::grid.grab()
}

mutid.v <- unname(rep(x=mut.id.v, each=length(sigEpLim.id.v)))
sigEpLimid.v <- rep(x=sigEpLim.id.v, times=length(mut.id.v))
len <- unique(length(mutid.v), length(sigEpLimid.v))

col.v <- c(`-1`="lightblue", `0`="gray80", `1`="darkred")
breaks.v <- c(-1, -0.5, 0.5, 1)

for(loc.id in loc.id.v){
  
  hm.lst <- sapply(X=1:len, simplify=F, FUN=function(i){
    
    mut.id <- mutid.v[i]
    sigEpLim.id <- sigEpLimid.v[i]
    
    mx <- DF[DF$mut.id==mut.id & DF$loc.id==loc.id & DF$sigEpLim.id==sigEpLim.id,
             c("ind", "SIG.id", "trendnum")]
    mx <- reshape2::acast(data=mx, formula=SIG.id~ind, value.var="trendnum")
    
    par(cex.main=0.5)
    heatmap.2(x=mx, Rowv=T, Colv=F, dendrogram="row", 
              main=paste0(out.id, "_", mut.id, "_", sigEpLim.id,  "\n", loc.id),
              scale="none", trace="none", na.rm=F, na.color="darkblue", 
              margins=c(9, 7), col=col.v, breaks=breaks.v, cexRow=0.2, cexCol=1, 
              rowsep=seq(1:nrow(mx)), colsep=seq(1:ncol(mx)), sepwidth=c(0.01, 0.01),
              key=T, distfun = function(x) dist(x, method = "euclidean"))
    
    # CHECK HOW TO CAPTURE COMPLEX HEATMAP AS GROB
    #grid.grabExpr(
    #  ComplexHeatmap::Heatmap(matrix=mx, col=col.v, na_col="black",
    #                        cluster_columns=F, cluster_rows=T, 
    #                        row_dend_width=unit(10,"mm"), 
    #                        row_split=10, #column_split=seq(1:ncol(mx)),
    #                        row_names_gp=gpar(fontsize=10)
    #                        )
    #)                
    
    grab_grob()  
    
  })
  
  grid::grid.newpage()
  pdf(file=paste0(out.dir, "/", out.id, "_", loc.id, "_heatmap.pdf"), height=60, width=40)
  gridExtra::grid.arrange(grobs=hm.lst, ncol=length(sigEpLim.id.v), clip=TRUE)
  dev.off()
  
  rm(hm.lst, loc.id)
  
}

# rm(list=ls()); gc()