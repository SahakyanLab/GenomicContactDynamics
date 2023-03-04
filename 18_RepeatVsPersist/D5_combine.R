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
src.dir = out.dir = paste0(wk.dir, "/out_combine")
ijcount.file = paste0(wk.dir, "/min2Mb_ij_2Mb_Cp1To21_ijcount.txt")
### OTHER SETTINGS #############################################################
GROUP.CLASS <- list(
  not.transposon = c("Low_complexity", "RNA", "rRNA", "Satellite", "scRNA", 
                     "Simple_repeat", "snRNA", "srpRNA", "tRNA"),
  DNA.transposon = c("DNA", "DNA?"), 
  retro.transposon = c("LINE", "SINE", "LTR", "RC", "LINE?", "SINE?", "LTR?"),
  not.classified = c("Other", "Unknown", "Unknown?")
)
group.cols = setNames(object=c("#7E6148FF", "#3C5488FF", "#00A087FF", "black"), 
                      nm=names(GROUP.CLASS))
Cps = 1:21
Cps.forCalc = 19:21
src.id = paste0("Cp", Cps.forCalc[[1]], "To", tail(Cps.forCalc, n=1))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(cowplot)
library(data.table)
library(ggsci)
library(ggrepel)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Final data for plotting
COMB <- list()

# Choose elements to focus on based on fraction of persistent contacts with at 
# least 1 shared site

load(paste0(src.dir, "/minrep_metric_extremes_", src.id, ".RData"))
extremes <- df; rm(df)
setnames(extremes, old="values", new="minrep.nonzero.fr")

extremes <- extremes[!duplicated(extremes$repName),] # Two repNames classified in two fams and classes
COMB[["minrep.nonzero.fr"]] <- extremes

elm.order <- COMB[["minrep.nonzero.fr"]]$repName

# Append metric estimates

metrics <- c("min", "skew", "sum")
metric.len <- length(metrics)
for(m in 1:metric.len){
  
  metric.id <- paste0(metrics[[m]], "rep")
  load(paste0(src.dir, "/", metric.id, "_estimate_", src.id, ".RData"))
  df <- df[!duplicated(df$repName),]
  
  if( identical(elm.order, df$repName) ){
    
    estimate <- df[,"values", drop=F]
    rm(df)
    setnames(estimate, old="values", paste0(metric.id, ".est"))
    COMB[[metric.id]] <- estimate

  } else {
    stop("Estimate: Wrong element order.")
  }
  
}

# 

load(paste0(src.dir, "/hg19_site_coverage_length.RData"))
if( identical(elm.order, dimnames(mx)[[1]]) ){
  
  sitecov <- mx
  rm(mx)
  COMB[["sitecov"]] <- sitecov[,c("total.bp", "mean.len.bp", "copynum")]
  COMB$sitecov[,"copynum"] <- COMB$sitecov[,"copynum"] / sum(COMB$sitecov[,"copynum"]) 

} else {
  stop("Sitecov: Wrong element order.")
}

#

load(paste0(src.dir, "/minrep_contact_coverage.RData"))
if( identical(elm.order, dimnames(mx)[[1]]) ){
  
  ijcov <- mx
  rm(mx)
  COMB[["ijcov"]] <- ijcov[, as.character(Cps.forCalc), drop=F]
  COMB[["ijcov"]] <- data.frame(ijcov=rowSums(COMB[["ijcov"]], na.rm=T))
  
  ij.count <- setNames(as.numeric(readLines(ijcount.file)), nm=Cps)
  ij.count <- sum(ij.count[Cps.forCalc])
  COMB[["ijcov"]] <- COMB[["ijcov"]] / ij.count

} else {
  stop("ijcov: Wrong element order.")
}


#

names(COMB) <- NULL
COMB <- do.call("cbind", COMB)

## Filter

COMB <- COMB[COMB$plot.group != "not.transposon",] # 968 of 1395
COMB <- na.omit(COMB) # 937 off 968

filterElm = T
if(filterElm){
  
  #> unlist(lapply(COMB, FUN=function(x) any(is.na(x))))
  #repClass           repName minrep.nonzero.fr         repFamily        plot.group        minrep.est 
  #FALSE             FALSE              TRUE             FALSE             FALSE              TRUE 
  #skewrep.est        sumrep.est          total.bp       mean.len.bp           copynum             ijcov 
  #TRUE              TRUE             FALSE             FALSE             FALSE             FALSE 
  
  COMB <- COMB[COMB$minrep.nonzero.fr > 0.5, ] # 31 of 830
  
}

## PLOTS

COMB <- COMB[order(COMB$minrep.nonzero.fr, decreasing=F),]
COMB$repName <- factor(COMB$repName, levels=COMB$repName)

P.LST <- sapply(c("copynum", "minrep.nonzero.fr",
                  "sumrep.est", "minrep.est", "skewrep.est",
                  "ijcov", "mean.len.bp"), simplify=F, FUN=function(val.nme){
  
  df <- cbind.data.frame(repName=COMB$repName, 
                         values=COMB[[val.nme]],
                         plot.group=COMB$plot.group)
              
  p <- ggplot(data=df, aes(x=repName, y=values)) +
    #geom_col(aes(col=plot.group, fill=plot.group)) +
    geom_segment(aes(x=repName, xend=repName, y=0, yend=values), col="gray50", alpha=0.1) +
    geom_point(aes(col=plot.group), shape=1, size=2) + 
    geom_text_repel(aes(label=repName), size=0.5, box.padding=0.5) + 
    #geom_hline(yintercept=0.5, linetype="dashed", color="black", size=0.5) + 
    scale_color_manual(values=group.cols) + 
    labs(y=val.nme) + 
    bgr1 + 
    theme(axis.title.x=element_blank(), axis.text.x=element_blank(), axis.ticks.x=element_blank(),
          axis.line.x=element_line(color="black"), 
          axis.title.y=element_text(size=5), axis.line.y=element_blank(),
          panel.background=element_blank(), legend.position="none",
          strip.background=element_blank(), strip.text=element_blank()) +
    facet_grid(.~plot.group) 
  
  return(p)
  
})

row.num = length(P.LST)
col.num = length(P.LST) / row.num
p.arr <- plot_grid(plotlist=P.LST, nrow=row.num, ncol=col.num, align="hv", byrow=T)
save_plot(filename=paste0(out.dir, "/out_filterElm", filterElm, ".pdf"), plot=p.arr,
          base_width=col.num * 5, base_height=row.num * 1.5)

writeLines(as.character(COMB$repName), paste0(out.dir, "/out_filterElm", filterElm, "_repName.txt"))

# rm(list=ls()); gc()