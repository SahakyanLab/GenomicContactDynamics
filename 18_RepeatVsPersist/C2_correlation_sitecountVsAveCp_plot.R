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
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")

# ELMTISSDYN Only for getting subfamily names
elm.dir = paste0(wk.dir, "/out_HicRepeatHeatmap/subfam_minrep_atleast2sumrep")
cor.dir = paste0(wk.dir, "/out_sitesVsmeanCp_binsWithNoSiteIncluded") #binsWithSiteOnlyIncluded")
out.dir = paste0(wk.dir, "/out_correlation_sitecountVsAveCp_plot/binsWithNoSiteIncluded") #binsWithSiteOnlyIncluded")
repTranspoType.file = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/17_RepeatAge/out_cleanAgeRank/agerank_summary_sheet1.csv")
### OTHER SETTINGS #############################################################
elm.id = "chrALL_min2Mb_ElmTissDyn_GiorPubl"
cor.id = "wmeanCp" # "raw" | "norm" | "fc"
alpha.val = 0.05
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
library(data.table)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.id <- paste0(elm.id, "_", cor.id, "_alpha", alpha.val)

load(paste0(elm.dir, "/", elm.id, ".RData"))
element <- rownames(ELMTISSDYN[[1]])
element.len <- length(element)
rm(ELMTISSDYN)

# Extract correlation data

mx <- matrix(data=NA, ncol=6, nrow=element.len, 
             dimnames=list(element, c("pear.coef", "pear.pval", "spea.coef", "spea.pval",
                                      "kend.coef", "kend.pval")))
for(elm.ind in 1:element.len){
  
  elm <- element[[elm.ind]]
  cor.path <- paste0(cor.dir, "/", cor.id, "_", elm, "_cortest.RData_cortest.RData")
  if( file.exists(cor.path) ){
    load(cor.path)
    
    cor.types <- setdiff(names(TEST), "alt")
    for(typ in cor.types){
      
      mx[elm.ind, paste0(typ, ".coef") ] <- unname( TEST[[typ]][["estimate"]] )
      mx[elm.ind, paste0(typ, ".pval") ]  <- unname( TEST[[typ]][["p.value"]] )
      
    }
    
    message(paste0(elm.ind, ": Correlation data obtained."))
    
  } else {
    message(paste0(elm.ind, ": No correlation data."))
  }
  
}

#

df <- as.data.frame(mx)
rm(mx)
df$element <- rownames(df)
df$Rank <- 1:element.len

# # Add clustering result
# clusters <- readLines(con=paste0(clust.dir, "/", clust.id, ".txt"))
# 
# is_nmes <- grepl(x=clusters, pattern=">Cluster", fixed=T)
# clusters.lst <- lapply(X=clusters[!is_nmes], FUN=function(str) strsplit(x=str, split=";")[[1]])
# 
# clust.nmes <- clusters[is_nmes]
# clust.nmes <- gsub(clust.nmes, pattern=">", replacement="")
# names(clusters.lst) <- clust.nmes
# clust.df <- stack(clusters.lst)
# colnames(clust.df) <- c("element", "cluster.name")
# 
# clust.df$cluster.name <- as.character(clust.df$cluster.name)
# indInclustdf <- match(x=df$element, table=clust.df$element)
# 
# df$cluster <- NA
# df$cluster <- clust.df$cluster.name[indInclustdf]
# 
# rm(clusters, is_nmes, clusters.lst, clust.nmes, clust.df, indInclustdf)

# Add if retrotransposon or DNA transposon

repTranspoType.df <- read.csv(file=repTranspoType.file, header=T)
if( identical(repTranspoType.df$repName, df$element) ){
  df <- cbind(df, repTranspoType.df[,"repTranspoType", drop=F])
}

# Plot

fixed.cols <- c("repTranspoType", "Rank") #c("cluster", "repTranspoType", "Rank")
for(typ in cor.types){
  
  df.tmp <- df[ ,c(paste0(typ, ".coef"), paste0(typ, ".pval"), fixed.cols) ]
  data.table::setnames(df.tmp, old=paste0(typ, ".coef"), "value")
  data.table::setnames(df.tmp, old=paste0(typ, ".pval"), "pval")
  
  # P-value filtering
  df.tmp$sig <- 0
  is_sig <- df.tmp$pval < alpha.val
  df.tmp$sig[is_sig] <- 1
 
  #
  
  plot.title <- paste0(out.id, "_", typ, "_pointswithgrayshadowarenotsig>=alpha", alpha.val)
  out.name <- paste0(out.id, "_", typ)
  
  p <- ggplot(data=df.tmp, aes(x=Rank, y=value)) +
    geom_hline(yintercept=0, lty="dashed", alpha=0.5) + 
    geom_point(data=df.tmp[df.tmp$sig==0,], size=10, shape=20, colour="gray70") + 
    geom_point(size=3.5, alpha=0.7, aes(shape=repTranspoType)) +
    scale_x_reverse() + 
    scale_y_continuous(limits=c(-0.3, 0.3)) +
    scale_shape_manual(values=c(4,19)) + 
    #scale_colour_manual(values=c("#541352FF", "#2f9aa0FF", "#ffcf20FF")) + 
    labs(x="Rank", y="Cor", title=plot.title) + 
    bgr2 + 
    theme(legend.position="bottom",legend.text=element_text(size=5),
          legend.title=element_blank(), axis.title.x=element_blank(),
          plot.title=element_text(size=5), aspect.ratio=1.5) +
    coord_flip()
  
  ggsave(filename=paste0(out.dir, "/", out.name, "_corplot.pdf"),
         width=10, height=10, unit="in", plot=p)
  
  save(df.tmp, file=paste0(out.dir, "/", out.name, "_corplot.RData"))
  
  rm(df.tmp, is_sig, plot.title)
  
  message(paste0(typ, ": Plotted!"))
  
}

# rm(list=ls()); gc()
