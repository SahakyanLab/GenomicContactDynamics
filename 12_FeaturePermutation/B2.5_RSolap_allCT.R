################################################################################
# Consolidate RSolap per tissue into one plot.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(wk.dir, "/out_RSolap_perCT/CpAllCsmatch_HiCNorm")
out.dir = paste0(wk.dir, "/out_RSolap_allCT/CpAllCsmatch_HiCNorm")
### OTHER SETTINGS #############################################################
data.id = "min2Mb_allCT" 

# UpSet plot
ylim.v = c(0,1)
at.v = seq(from=0.2, to=0.8, by=0.2)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ComplexHeatmap)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(data.id, "..."), quote=FALSE)
#-------------------Extract info, including boxplot data, from UPSMX
load(file=paste0(data.dir, "/", data.id, "_RS_upsmx.RData"))
ct.v <- setdiff(rownames(UPSMX$CS), "allCT")
RSid.v <- colnames(UPSMX$SS)
comb.nme <- colnames(UPSMX$CS)[ colSums(UPSMX$CS[ct.v,])>0 ]
totreg <- UPSMX$totunionreg

UPSMX$CS <- sapply(X=comb.nme, simplify=FALSE, FUN=function(nme){
  return( unname(UPSMX$CS[ct.v,nme]) )
})
UPSMX$SS <- sapply(X=RSid.v, simplify=FALSE, FUN=function(RSid){
  return( unname(UPSMX$SS[ct.v,RSid]) )
})
#-------------------UpSet plot
load(file=paste0(data.dir, "/", data.id, "_RS_upsetplot.RData"))
UPSOBJ <- UPSOBJ[comb_name(UPSOBJ)%in%comb.nme]
cs <- comb_size(UPSOBJ)/totreg
ss <- set_size(UPSOBJ)/totreg
comb.od <- order(-comb_degree(UPSOBJ))

pdf(file=paste0(out.dir, "/", data.id, "_upsetplot_comb.pdf"), width=10, height=10)
ht <- UpSet(UPSOBJ, 
            set_order=match(x=set_name(UPSOBJ), table=RSid.v), 
            comb_order=comb.od,
            top_annotation=HeatmapAnnotation(height=unit(2, "cm"), annotation_name_side="left",
              cs=anno_barplot(cs, ylim=ylim.v, border=FALSE, gp=gpar(fill="black"),
                              axis_param=list(at=at.v))
            ),
            bottom_annotation=HeatmapAnnotation(height=unit(2, "cm"), annotation_name_side="left",
              f=anno_boxplot(UPSMX$CS[comb_name(UPSOBJ)], outline=TRUE, 
                             ylim=ylim.v, axis_param=list(at=at.v), border=TRUE)
            ),
            right_annotation=rowAnnotation(width=unit(2, "cm"),
              ss=anno_barplot(ss, baseline=0, border=FALSE, gp=gpar(fill="black"),
                              axis_param=list(at=at.v))
            ),
            left_annotation=rowAnnotation(width=unit(2, "cm"),
              f=anno_boxplot(UPSMX$SS[RSid.v], outline=TRUE, ylim=ylim.v, 
                             axis_param=list(at=at.v), border=TRUE)
            ),
            row_title=paste0(data.id, "\ndistinctmode\ntotalunionregion=", 
                             totreg), show_row_names=TRUE
)
ht <- draw(ht)
c.od <- column_order(ht)
r.od <- row_order(ht)
decorate_annotation("cs", {
  cs.lab <- cs[c.od]; cs0.TF <- cs.lab==0
  cs.lab <-  format(cs.lab, digits=2) 
  cs.lab[cs0.TF] <- "0"
   grid.text(label=cs.lab, x=1:length(cs), 
             y=unit(cs[c.od], "native") + unit(2, "pt"), 
            default.units="native", just=c("center", "bottom"), 
            gp=gpar(fontsize=5, col="black"), rot=0)
})
decorate_annotation("ss", {
  grid.text(label=format(ss[r.od], digits=3, scientific=FALSE), 
            y=length(ss):1, x=unit(ss[r.od], "native"), 
            default.units="native", just=c("center", "bottom"), 
            gp=gpar(fontsize=5, col="white"), rot=90)
})
dev.off()

# rm(list=ls()); gc()



