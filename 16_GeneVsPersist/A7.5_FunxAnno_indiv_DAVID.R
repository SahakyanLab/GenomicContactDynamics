################################################################################
# Plot DAVID functional annotation clusters of terms enriched in given gene set.
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
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}

david.dir = paste0(wk.dir, "/out_DAVID/funxAnnoClust_sampleGeneList")
out.dir = paste0(wk.dir, "/out_FunxAnno_indiv_DAVID")
### OTHER SETTINGS #############################################################
data.id = "min2Mb_LTr_ALL_cp_21_name2"
seed.v = c(287, 587, 754)
N = 2999
stringency = "medium"

# Get only top 10 clusters with highest enrichment score
Nclust = 5

putlabel = T
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(data.id, "_top", Nclust, "highestEASEclusters_seed",
                 paste(seed.v, collapse="_"), "_", N, "sampled_ftc_", stringency,
                 "_putlabel", putlabel)
  
FTC <- list()
ANME <- list() # Annotation cluster name and enrichment score
for( seed in as.character(seed.v) ){
  
  FTC[[seed]] <- readLines(con=paste0(david.dir, "/", data.id, "_seed", seed, 
                                      "_", N, "sampled_ftc_", stringency, ".txt"))
  
  # Remove empty strings usually before and after cluster name
  FTC[[seed]] <- FTC[[seed]][FTC[[seed]]!=""]
  
  # Remove headers ("Category\tTerm\tCount") per cluster
  headr.TF <- grepl(x=FTC[[seed]], pattern="Category\tTerm\tCount", fixed=T)
  headr <- FTC[[seed]][headr.TF][1] # Store
  headr <- strsplit(x=headr, split="\t", fixed=T)[[1]]
  FTC[[seed]] <- FTC[[seed]][!headr.TF]
  
  # Indices of cluster names
  cind.v <- which(grepl(x=FTC[[seed]], pattern="Annotation Cluster", fixed=T))
  
  # Check that number of headers same as number of clusters
  if( sum(headr.TF)!=length(cind.v) ){
    stop(paste0(seed, ": Checkpoint 1."))
  } 
  
  # Get top 10 clusters
  FTC[[seed]] <- FTC[[seed]][ 1:(cind.v[Nclust+1]-1) ]
  
  # Get cluster name, enrichment score and number of member terms
  anme <- paste(FTC[[seed]][ cind.v[1:Nclust] ], 
                # Number of terms per cluster
                ( diff(cind.v[1:(Nclust+1)]) - 1 ), sep="\t")
  anme <- gsub(x=anme, pattern="Annotation Cluster", replacement="", fixed=T)
  anme <- gsub(x=anme, pattern="Enrichment Score:", replacement="", fixed=T)
  
  # Transform FTC to dataframe
  FTC[[seed]] <- FTC[[seed]][ -cind.v[1:10] ]
  FTC[[seed]] <- read.delim(text=FTC[[seed]], sep="\t", stringsAsFactors=F, header=F)
  colnames(FTC[[seed]]) <- headr
  
  # Add cluster info to dataframe
  anme <- read.delim(text=anme, sep="\t", stringsAsFactors=F, header=F)
  colnames(anme) <- c("Cluster", "EASE", "Nterms")
  if( sum(anme$Nterms)==length(FTC[[seed]][[1]]) ){
    
    for( nme in c("Cluster", "EASE") ){
      
      anme[[nme]] <- stringr::str_trim(string=anme[[nme]], side="both")
      anme[[nme]] <- as.numeric(as.character(anme[[nme]]))
      FTC[[seed]][[nme]] <- unlist(mapply(FUN=rep, x=anme[[nme]], each=anme$Nterms))
      
    }
    
  } else {
    stop(paste0(seed, ": Checkpoint 2."))
  }
  
  FTC[[seed]]$seed <- as.character(seed)
  
  rm(anme, cind.v)

}

FTC <- do.call("rbind", FTC)
save(FTC, file=paste0(out.dir, "/", out.id, "_DAVIDplot.RData"))

# Take only terms significant in at least one sample;
# Use Benjamini to match the BH p-adjusted values plotted for KEGG and GO
#FTC <- FTC[ FTC$Term%in%FTC$Term[FTC$Benjamini<0.05], ]
FTC <- FTC[order(FTC$Benjamini, decreasing=F), ]

FTC$Term <- factor(x=FTC$Term, levels=rev(unique(FTC$Term)))
FTC$Cluster <- factor(x=as.character(FTC$Cluster), 
                      levels=as.character(sort( unique(FTC$Cluster)), decreasing=F ))


clust.col <- colorRampPalette(RColorBrewer::brewer.pal(Nclust, "Spectral"))(Nclust)
names(clust.col) <- levels(FTC$Cluster)
ptitle <- Xlab <- NULL
if(putlabel){
  
  ptitle <- paste0(out.id, "_Cluster", Clust, "_reddashedlineis-log10(0.05)")
  Xlab <- bquote(bold( "-log10("~"p-value"^"adj Benjamini"~")" ))
  
}

p.lst <- list()
for( Clust in levels(FTC$Cluster) ){
  
  Clust <- as.character(Clust)
  p.lst[[Clust]] <- ggplot(data=FTC[as.character(FTC$Cluster)==Clust,], 
                           aes(x=-log10(Benjamini), y=Term) ) +
    geom_point(fill=clust.col[[Clust]], aes(size=Count), colour="black", shape=21) + 
    geom_vline(xintercept=-log10(0.05), linetype="dashed", colour="red", size=0.7) +
    scale_x_continuous(limits=c(0,10)) + 
    guides(colour="legend") +
    labs(title=ptitle, size="Gene\ncount", x=Xlab, y=NULL) + 
    bgr5 + 
    theme(panel.background=element_rect(colour="gray22", size=1, fill=NA),
          plot.title=element_text(size=5),
          axis.text.y=element_text(size=15),
          axis.text.x=element_blank(),
          legend.title=element_text(size=20), 
          legend.text=element_text(size=20)) +
    facet_grid(.~seed)  
  
  if(!putlabel){
    
    p.lst[[Clust]] <- p.lst[[Clust]] + 
      theme(strip.text=element_blank())
    
  }
  
  if( Clust==as.character(Nclust) ){
    
    p.lst[[Clust]] <- p.lst[[Clust]] + 
      theme(axis.text.x=element_text(size=20, face="bold"))

  }
  
}

seed.v.len <- length(seed.v)
p.arr <- ggpubr::ggarrange(plotlist=p.lst, nrow=Nclust, ncol=1, align="hv")
ggpubr::ggexport(p.arr, width=20, height=5*Nclust,
                 filename=paste0(out.dir, "/", out.id, "_DAVIDplot.pdf"))

# rm(list=ls()); gc()




