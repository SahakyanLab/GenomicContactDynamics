################################################################################
# Plot per location, per filter, combining all signatures
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
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
plot.dir = paste0(wk.dir, "/out_plotsVsCp")
out.dir = paste0(wk.dir, "/out_signatureCombinedPlot")
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
library(yarrr)
source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
makePlot <- function(df, calc, yint=NULL, plot.id){
  
  yarrr::piratepal("basel")
  
  if( !is.numeric(df$ind) ){
    stop("Plot making: ind not numeric as expected.")
  }
 
  df$pval <- c(`0`="ns", `1`="sig")[ as.character(as.numeric(df$pval<0.05)) ]
  df$pval <- factor(df$pval, levels=c("ns", "sig"))
  shape.v <- c(8, 1)
  names(shape.v) <- levels(df$pval)
  
  col.v <- colorRampPalette(yarrr::piratepal("basel"))(length(unique(df$SIG.id)))
  df$SIG.id <- factor(df$SIG.id, levels=unique(c("nosampfilter", sort(df$SIG.id))))
  # Assign colour black to nosampfilter
  col.v[1] <- "black"
    
  p <- ggplot(data=df, aes(x=ind, y=values)) +
    geom_line(aes(colour=SIG.id), alpha=0.3) +
    geom_point(aes(colour=SIG.id, shape=pval), size=2) +
    scale_x_continuous(breaks=unique(sort(df$ind))) + 
    scale_colour_manual(values=col.v) +
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
  
  pdata.id <- paste0(gcb, "_", data.id, "_", src.id, "_", mut.id, "_", calc, "_",
                     aggregatefunx)
  load(file=paste0(plot.dir, "/", pdata.id, "_scatplot.RData"))    # P.DF
  load(file=paste0(plot.dir, "/", pdata.id, "_FC_scatplot.RData")) # PFC.DF
  
  DF[[paste0(mut.id, calc)]] <- cbind.data.frame(PFC.DF, values.abs=P.DF$values,
                                                 stringsAsFactors=F)
  
  rm(P.DF, PFC.DF); gc()
  
} # mut.id.v for loop end

DF <- do.call("rbind.data.frame", c(DF, stringsAsFactors=F))
rownames(DF) <- NULL

# Separate nosampfilter rows
nosamp.df <- DF[DF$sigEpLim.id=="nosampfilter",]
nosamp.df <- nosamp.df[!duplicated(nosamp.df[,-2]),]
nosamp.df$SIG.id <- "nosampfilter"

DF <- DF[DF$sigEpLim.id!="nosampfilter",]

loc.id.v <- unique(DF$loc.id)
sigEpLim.id.v <- unique(DF$sigEpLim.id)
sigEpLim.id.v.len <- length(sigEpLim.id.v)
mut.id.v <- unique(DF$mut.id)
mut.id.v.len <- length(mut.id.v)

s.v <- rep(x=sigEpLim.id.v, times=length(mut.id.v))
s.v <- s.v[s.v!="nosampfilter"]
s.v.len <- length(unique(s.v))
m.v <- rep(x=mut.id.v, each=length(sigEpLim.id.v[sigEpLim.id.v!="nosampfilter"]))
m.v.len <- length(unique(m.v))
len <- unique(length(s.v), length(m.v))

for(loc.id in loc.id.v){
  
  p.lst <- list()
  for(i in 1:len){
    
    mut.id <- m.v[i]
    sigEpLim.id <- s.v[i]
    
    df <- DF[DF$loc.id==loc.id & DF$sigEpLim.id==sigEpLim.id & DF$mut.id==mut.id,]
    df <- rbind.data.frame(df, nosamp.df[nosamp.df$loc.id==loc.id & nosamp.df$mut.id==mut.id,])
    p.lst[[paste0(loc.id, "_", i)]] <- makePlot(df=df, 
                                                calc=calc, yint=0,
                                                plot.id=paste0(out.id, "\n", loc.id, "_", sigEpLim.id, "_", mut.id))
    rm(df)
    
  }
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=m.v.len, ncol=s.v.len, legend=NULL)
  ggexport(p.arr, height=10*m.v.len, width=10*s.v.len, 
           filename=paste0(out.dir, "/", out.id, "_", loc.id, ".pdf"))
  rm(p.lst, mut.id, sigEpLim.id)
  
}

# rm(list=ls()); gc()