################################################################################
# Characterize the repeat clusters from the RepeatsVsPersist analyses 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/18_RepeatAge"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
repFilePath = paste0(wk.dir, "/out_addToSummary/hg19repeats_repName.RData")
out.dir = paste0(wk.dir, "/out_clusters")
### OTHER SETTINGS #############################################################
rank.v = c("Giordano364rank", "GiorPubl372rank", "Publrank")
plotOnly = FALSE
agerankPlot = FALSE
CNplot = FALSE
RepTypeOriginPlot = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggsci)
library(viridis)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Load REPEAT.MX (data.frame)
load(file=repFilePath)
REPEAT.MX$copyNumber <- as.numeric(REPEAT.MX$copyNumber)
# Count repName/subfamily per cluster
repSFperClust <- table(REPEAT.MX[,"cluster"])
# Percentage of sites per cluster relative to total repeat sites
percSitesperClust <- by(data=REPEAT.MX$copyNumber, INDICES=REPEAT.MX$cluster,
                        FUN=sum)
percSitesperClust <- percSitesperClust/sum(REPEAT.MX$copyNumber)*100
percSitesperClust <- format(x=percSitesperClust, digits=2, scientific=TRUE)
percSitesperClust <- paste(percSitesperClust, "%", sep="")
names(percSitesperClust) <- names(repSFperClust)
#-----------------------------
# Boxplot of repeat subfamily copynumbers per cluster
if(CNplot==TRUE){
  
  df <- REPEAT.MX[,c("cluster", "copyNumber")]
  #-------------------
  # Plot with outliers
  p <- ggplot(data=df, aes(x=cluster, y=copyNumber)) +
    stat_boxplot(geom="errorbar", width=0.4) + 
    geom_boxplot(fill="#55bde6") +
    annotate("text", x=names(repSFperClust), y=-6000, 
             label=repSFperClust, size=3) +  
    annotate("text", x=names(percSitesperClust), y=-10000, 
             label=percSitesperClust, size=3) +
    labs(x=NULL,
         y="Copy Number",
         title="hg19repeats_Clusters_copyNumber") +
    scale_x_discrete(breaks=names(repSFperClust)) +
    bgr1 +
    theme(axis.text.x=element_text(size=15, angle=360, colour="black"))
  ggsave(file=paste0(out.dir, "/plot_hg19_clustersVsCopyNumber_outl.pdf"), 
         units="in", width=10, height=10, plot=p)
  #-------------------
  # Zoom in to plot
  # Compute the extreme upper and lower whiskers to scale plot after 
  ylim.val <- boxplot.stats(x=df$copyNumber)$stats[c(1, 5)]
  # Make new plot
  p <- p + 
    coord_cartesian(ylim = c(-1000, ylim.val[2]*10)) +
    annotate("text", x=names(repSFperClust), y=-1000, 
             label=repSFperClust, size=3) +  
    annotate("text", x=names(percSitesperClust), y=-2000, 
             label=percSitesperClust, size=3) 
  ggsave(file=paste0(out.dir, "/plot_hg19_clustersVsCopyNumber_zoomed.pdf"), 
         units="in", width=10, height=10, plot=p)
  rm(df, p, ylim.val)
  
}
#-----------------------------
# Plot ClustersVsRepType and ClustersVsAnimalorigin
if(RepTypeOriginPlot==TRUE){
  
  # Animalorigin: Warning message:
  # Removed 301 rows containing non-finite values (stat_sum).
  # Cause: no cluster contain Euteleostemi
  # No NAs for repType
  for(col in c("repType", "origin", "repClass")){
    
    if(plotOnly==FALSE){
      df <- REPEAT.MX[,c("cluster", col, "copyNumber")]
      df[[col]][ is.na(df[[col]]) ] <- "Unknown"
      df <- aggregate(x=df$copyNumber, by=list(df$cluster, df[[col]]), FUN=sum)
      colnames(df) <- c("cluster", "group", "copyNumber")
      save(df, file=paste0(out.dir, "/plot_hg19_clustersVs", col, ".RData"))
    } else {
      load(file=paste0(out.dir, "/plot_hg19_clustersVs", col, ".RData"))
    }

    clust.CN <- aggregate(x=df$copyNumber, by=list(df$cluster), FUN=sum)
    clust.CN <- paste(paste(clust.CN$Group.1, clust.CN$x, sep=":"), collapse="_")
    coul <- colorRampPalette(ggsci::pal_npg("nrc")(9))(length(unique(df$group)))
    
    p <- ggplot(data=df, aes(x=cluster, y=copyNumber, fill=group)) +
      geom_col(position="fill") +
      scale_fill_manual(values=coul) +
      labs(title=paste0("hg19repeats_Clusters_", col, "\nCNpergroup_", clust.CN), 
           x=NULL, y="Fraction of insertion sites", 
           fill=col) + 
      bgr2 +
      theme(axis.text.x=element_text(angle=90, hjust=1))
   
    ggsave(file=paste0(out.dir, "/plot_hg19_clustersVs", col, ".pdf"),
           plot=p, units="in", width=10, height=10)
    
  }
  
}
#-----------------------------
# plot clustersVsAgeranks
if(agerankPlot==TRUE){
  
  rank.v.len <- length(rank.v)
  for(k in 1:rank.v.len){
    
    if(plotOnly==FALSE){
      df <- REPEAT.MX[REPEAT.MX[[ rank.v[k] ]]!=0,
                      c("copyNumber", "cluster", rank.v[k])]
      melt1 <- do.call( "c", mapply(FUN=rep, df[,"cluster"], 
                                    df[,"copyNumber"]))
      melt2 <- do.call( "c", mapply(FUN=rep, df[,rank.v[k]], 
                                    df[,"copyNumber"]))
      df <- data.frame(cluster=melt1, agerank=melt2, stringsAsFactors=FALSE)
      save(df, file=paste0(out.dir, "/plot_hg19_clustersVs", 
                           rank.v[k], ".RData"))
    } else {
      load(file=paste0(out.dir, "/plot_hg19_clustersVs", rank.v[k], ".RData"))
    }
    numTEs <- length(unique(df$agerank))
    p <- ggplot(data=df, aes(x=agerank)) +
      # y=..count.. to reflect the difference in number of sites each cluster
      #geom_density( alpha=0.4, aes(fill=factor(cluster), y=..count..) ) +
      geom_density( alpha=0.4, aes(fill=factor(cluster)) ) +
      guides(fill=FALSE) + 
      scale_x_continuous(limits=c(1, numTEs), breaks=c(0,100,200,249,344,364)) + 
      scale_fill_manual(values=viridis(n=3)[1:2]) + 
      labs(title=paste0("hg19repeats_", rank.v[k]),
           x=paste0(numTEs, " TEs in chronological order"),
           y="Density") + 
      facet_grid(cluster ~ .) + 
      bgr1 + 
      theme(strip.text.y=element_text(size=12, angle=90, face="bold"),
            aspect.ratio=0.2) 
    ggsave(file=paste0(out.dir, "/plot_hg19_clustersVs", rank.v[k], ".pdf"),
           plot=p, units="in", width=10, height=10)
    rm(df, p)
    
  }
  
}

# rm(list=ls()); gc()

