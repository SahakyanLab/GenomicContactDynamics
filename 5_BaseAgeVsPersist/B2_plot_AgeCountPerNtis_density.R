################################################################################
#Density plot to show distribution of the number of bases per bin with age data

################################################################################
#functions to source
#lib = "/Users/ltamon/ProjectBoard/lib"
lib = "/t1-home/icbml/ltamon/ProjectBoard/lib"
#main directory for the task
#objective.dir = "/Users/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/5_BaseAgeVsPersist"
data.dir = paste0(objective.dir, "/out_AgeCountPerNtis")  
output.dir = paste0(objective.dir, "/out_plot_AgeCountPerNtis_density")

library(foreach)
#library(ggpubr) #also loads ggplot2
library(ggplot2)

#2(2MB gap) or "05"(0.5 MB gap), refers to minimum gap accepted to classify 
#a contact, two points should be far enough to filter for contacts within a TAD
gcb.v = c("2", "05")

chr.v = c(1:22, "X")

################################################################################
source(paste0(lib, "/ggLayers_persist.R"))
#source(paste0(lib, "/multiplot.R"))

densplot <- function(df=df){
  ggplot(data=as.data.frame(df), aes(count)) +
    geom_density(position = "stack", aes(y=..count.., fill=factor(ntis))) +
    geom_vline(aes(xintercept=mean(count)),
               color="black", linetype="dashed", size=0.5) +
    labs(fill="Persistence score") +
    scale_fill_manual(values=PersistScoreColour(ntis.uniq) ) +
    ggtitle(paste0("chr", chr, " _min", gcb, "Mb")) +
    xlab("Number of base age data per bin") +
    #include number of datapoints per ntis
    ylab(paste0("Density (count)")) +
    theme.persist #from sourced "ggLayers_persist.R"
}

################################################################################
foreach(gcb=gcb.v, .inorder=TRUE) %do% {
  myplots.list <- list()
  for(chr in chr.v){
      load(file=paste0(data.dir,"/plot_min", gcb, "Mb_chr", 
                       chr, "_AgeCountPerNtis.RData"))
      ntis.uniq <- sort( as.numeric( unique(AGECOUNT.NTIS.MX[,"ntis"]) ),
                         decreasing=FALSE)
    #N <- nrow(AGECOUNT.NTIS.MX)
    p <- densplot(df=AGECOUNT.NTIS.MX)
    #myplots.list[[paste0("chr", chr)]] <- p
    ggsave(file=paste0(output.dir, "/dens_each_min", gcb, "Mb_chr", chr, 
                       "_AgeCountPerNtis.jpeg"), plot=p)
  }
  #pdf(file=paste0(output.dir, "/plot_min", gcb, "Mb_chr", paste(chr.v, collapse=""), 
  #                "_AgeCountPerNtis.pdf"), width=13, height=8.5)
  #multiplot(plotlist=myplots.list, cols=3)
  #dev.off()
  #myplots.arranged <- ggarrange(plotlist=myplots.list, nrow=2, ncol=3,
  #                              common.legend=TRUE, legend="right")
  #ggexport(myplots.arranged, width=13, height=8.5,
  #         filename=paste0(output.dir, "/densplot_min", gcb, 
  #                         "Mb_chr", paste(chr.v, collapse=""),
  #                         "_AgeCountPerNtis.pdf"))
}

#rm(list=ls())

