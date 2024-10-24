################################################################################
# Plot showing the type (DNA/retrotransposon), family and cluster of repeat
# subfamilies in the age ranks.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/17_RepeatAge"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
repFilePath = paste0(wk.dir, "/z_ignore_git/out_addToSummary/hg19repeats_repName.RData")
### OTHER SETTINGS #############################################################
out.dir = paste0(wk.dir, "/out_characterizeAgeRank")
rank.v = c("Giordano364rank", "GiorPubl372rank", "Publrank")
plotOnly = FALSE
col.id = "CptopCP3" #"cluster"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(viridis)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
################################################################################
# Load REPEAT.MX 
load(file=repFilePath)
rank.v.len <- length(rank.v)

for(k in 1:rank.v.len){
  if(plotOnly==FALSE){
    df <- REPEAT.MX[REPEAT.MX[[ rank.v[k] ]]!=0, 
                    c(rank.v[k], "repName", "repFamily", "repType", col.id)]
    setnames(x=df, old=rank.v[k], new="rank")
    df$rank <- as.numeric(df$rank)
    df <- within(data=df, {
      repTranspoType <- NA
      repTranspoType[repType %in% c("non-LTR", "LTR")] <- "Retrotransposon"
      repTranspoType[repType=="DNAtransposon"] <- "DNAtransposon"
    })
    rownames(df) <- NULL

    df <- df[order(df$rank, decreasing=FALSE),]
    write.csv(df, quote=FALSE, row.names=FALSE,
              file=paste0(out.dir, "/plot_", rank.v[k], "repFamilies.csv"))
    save(df, file=paste0(out.dir, "/plot_", rank.v[k], "repFamilies.RData"))
  } else {
    load(file=paste0(out.dir, "/plot_", rank.v[k], "repFamilies.RData"))
  }
  
  df <- df[order(df$repFamily, decreasing=T),]
  
  ggplot(data=df, aes(x=rank, y=repFamily)) +
    geom_point( size=8, stroke=1, aes_string(colour=col.id,
                                             shape="repTranspoType") ) +
    guides(size="none") +
    labs(colour=NULL, shape=NULL, main=rank.v[k], y="Repeat Family", 
         x="TE subfamilies in chronological order") +
    scale_shape_manual(values=c(4,1)) +
    #scale_size_manual(values=c(4,4)) + 
    scale_y_discrete(limits=unique(df$repFamily)) + 
    scale_color_manual(values=viridis(n=3)[1:length(unique(REPEAT.MX[[col.id]]))]) +
    bgr2 + 
    theme(panel.grid.major.y=element_line(linetype="dashed",
                                          size=0.5, colour="grey80"),
          legend.text=element_text(size=25, face="bold"))
 
  ggsave(filename=paste0(out.dir, "/plot_", rank.v[k], "repFamilies_", col.id, ".pdf"),
         units="in", width=15, height=15)
  rm(df)
}

# rm(list=ls()); gc()
