################################################################################
# Plot number of sites family/subfamily per bin vs. mean Cp per bin
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")

out.dir = paste0(wk.dir, "/out_sitesVsmeanCp")
agerank.file = paste0(wk.dir, "/Repeat_rankingbyAge/repsubfam.csv")
### OTHER SETTINGS #############################################################
alpha = 0.05
#colorFam = c("Alu", "CR1", "ERV1", "ERVK", "ERVL", "ERVL-MaLR",
#             "L1", "L2", "MIR", "Other", "RTE")

#color.id = "retrotransposon" #"retrotransposon"

#colorFam = c("DNA", "hAT", "hAT-Blackjack", "hAT-Charlie", "hAT-Tip100", "MuDR",
#             "PiggyBac", "TcMar-Mariner", "TcMar-Tc2", "TcMar-Tigger")
            
#color.id = "dnatransposon" #"retrotransposon"

colorFam = NULL
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(viridis)
library(ggplot2)
library(ComplexHeatmap)
library(yarrr); base.pal <- yarrr::piratepal("basel")
source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
agerank <- read.csv(file=agerank.file, header=T, stringsAsFactors=F)
elm.v <- agerank[,"repName"]
elm.v.len <- length(elm.v)

for( wmCp in c("wmeanCp", "wmeanCp0") ){
  
  df <- matrix(data=NA, nrow=elm.v.len, ncol=5, 
               dimnames=list(elm.v, 
                             c("Rank", "Pear", "Pear.pval", "Spear", "Spear.pval"))) 
  df[,"Rank"] <- 1:elm.v.len
 
  for(elm.id in 1:elm.v.len){
    
    elm <- elm.v[elm.id]
    
    obj <- paste0(out.dir, "/", wmCp, "_", elm, "_cortest.RData_cortest.RData")
    
    if( file.exists(obj) ){
      
      load(obj)
      
      df[elm.id, -1] <- c(TEST$pear$estimate, TEST$pear$p.value,
                          TEST$spea$estimate, TEST$spea$p.value)
      
    } else {
      
      print(paste0(elm, ": No cor data to plot."))
      next
      
    }
  
  }
  
  out.name <- paste0("chrALL_min2Mb_GiorPubl372RankSiteCountPerBinVsBin", wmCp, "_cor")
  
  save(df, file=paste0(out.dir, "/", out.name, "_heatcorplot.RData"))
  
  df[df[,"Pear.pval"] > alpha, "Pear"] <- NA
  df[df[,"Spear.pval"] > alpha, "Spear"] <- NA
  
  plot.title <- paste0(out.name, "\n outOf", nrow(df), "elements-",
                       sum(is.na(df[,"Pear"])), "ZeroCorInPear_", sum(is.na(df[,"Spear"])), 
                       "ZeroCorInSpear_WithPval>", alpha, "SetTo0")
  
  # Heatmap
  
  pdf(file=paste0(out.dir, "/", out.name, "_corheatmap.pdf"), width=10, height=10)
  
  print(
    ComplexHeatmap::Heatmap(matrix=df[,c("Pear", "Spear")], col=viridis(n=299), na_col="gray50", 
                            cluster_columns=F, cluster_rows=F, 
                            row_names_gp=gpar(fontsize=2),
                            heatmap_legend_param=list(title=plot.title, at=seq(-0.3,0.3,length.out=7))
                            )
  )
  
  dev.off()
  
  # Scatter
  
  df <- as.data.frame(df)
  df$Pear.pval <- df$Spear.pval <- NULL
  df$repFamily <- agerank$repFamily
  #df$repFamily[ !df$repFamily %in%colorFam] <- NA
  
  if( !is.null(colorFam) ){
    df <- df[df$repFamily %in% colorFam,]
  } else {
    color.id <- "all"
  }
  
  plot.title <- paste0(out.name, "\n outOf", nrow(df), "elements-",
                       sum(is.na(df[,"Pear"])), "ZeroCorInPear_", sum(is.na(df[,"Spear"])), 
                       "ZeroCorInSpear_WithPval>", alpha, "SetTo0")
 
  df <- reshape2::melt(df, id.vars=c("repFamily", "Rank"))
  df$value[is.na(df$value)] <- 0
  
  p <- ggplot(data=df, aes(x=Rank, y=value)) +
    geom_hline(yintercept=0, lty="dashed", alpha=0.5) + 
    geom_point(data=df.tmp[df.tmp$sig==0,], size=10, shape=20, colour="gray70") + 
    geom_point(size=3.5, alpha=0.7, aes(colour=cluster, shape=repTranspoType)) +
    scale_x_reverse() + 
    scale_y_continuous(limits=c(-1, 1)) +
    scale_shape_manual(values=c(4,19)) + 
    scale_colour_manual(values=c("#541352FF", "#2f9aa0FF", "#ffcf20FF")) + 
    labs(x="Rank", y="Cor", title=plot.title) + 
    bgr2 + 
    theme(legend.position="bottom",legend.text=element_text(size=5),
          legend.title=element_blank(), axis.title.x=element_blank(),
          plot.title=element_text(size=5), aspect.ratio=1.5) +
    coord_flip()
  
  # ggplot(data=df, aes(x=Rank, y=value)) +
  #   geom_hline(yintercept=0) + 
  #   geom_point(size=2.5, shape=1, alpha=0.7, aes(colour=repFamily)) +
  #   #scale_y_continuous(limits=c(-0.3, 0.3)) + 
  #   scale_colour_manual(values=colorRampPalette(base.pal)(length(unique(df$repFamily)))) + 
  #   labs(x="Rank", y="Cor", title=plot.title) + 
  #   guides(colour=guide_legend(ncol=3)) +
  #   bgr2 + 
  #   facet_grid(.~variable) + 
  #   theme(legend.position="bottom",legend.text=element_text(size=15)
  #         )
  
  ggsave(filename=paste0(out.dir, "/", out.name, "_color", color.id, "_corplot.pdf"),
         width=10, height=10, unit="in", plot=p)
  
}

# rm(list=ls()); gc()