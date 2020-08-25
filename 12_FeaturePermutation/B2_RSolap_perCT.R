################################################################################
# Use UpSetR plot to visualise overlap of region sets amongst themselves.
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
#binmx.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc/binmx/out_bindata_1perc")
binmx.dir = paste0(wk.dir, "/binmx/out_bindata_Csmatch_HiCNorm")
out.dir = paste0(wk.dir, "/out_RSolap_perCT/CpAllCsmatch_HiCNorm")
### OTHER SETTINGS #############################################################
gcb = "min2Mb" 
chr.v = paste0("chr", c(1:22, "X"))
CT.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
         "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
# The three objects below should be in the same order
Cp.lst = list(21, 19:21, 1, 1:21, 1:21)
#Cs.lst = list(c(0.01,1), c(0.01,1), c(0.01,1), c(0.01,1), 0.01)
Cs.lst = list(c(2,1), c(2,1), c(2,1), c(2,1), 2)
RSid.v = c("Cp21", "CptopCP3", "Cp1", "CpAll", "CpAllCsmatch")

# UpSet plot axes
ylim.v = c(0, 80000)
at.v = seq(from=ylim.v[1], to=ylim.v[2], by=20000)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ComplexHeatmap)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(paste0(gcb, "..."), quote=FALSE)

CT.v <- sort(unique(CT.v))
Cp.lst <- lapply(X=Cp.lst, FUN=function(x){
  sort(unique(x))
})
Cs.lst <- lapply(X=Cs.lst, FUN=function(x){
  (unique(x))
})
names(Cp.lst) <- names(Cs.lst) <- RSid.v
#---------------------------------------
# Prepare contacting region coordinates filtering by Cp.v
if( is.null(chr.v) ){ chr.v <- paste0("chr", c(1:22, "X")) }
BIN.MX <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
  load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
  BIN.MX <- cbind.data.frame(chr=chr, BIN.MX, stringsAsFactors=FALSE)
})
BIN.MX <- do.call("rbind", BIN.MX)
BIN.MX <- data.matrix(BIN.MX[,!colnames(BIN.MX)%in%c("chr", "start", "end")])
BIN.MX  <- BIN.MX[rowSums(BIN.MX)>0,]
uni.set <- dimnames(BIN.MX)[[1]]; totreg <- length(uni.set)
if( any(duplicated(dimnames(BIN.MX))) ){ stop("Duplicated region names.") }
#---------------------------------------
CS <- SS <- list()
pdf(file=paste0(out.dir, "/", gcb, "_RS_upsetplot.pdf"), width=10, height=10)
for( ct in c(CT.v, "allCT") ){
  
  out.name <- paste0(gcb, "_", ct, "_RS_upsetplot")
  RS.lst <- list()
  for(RS in RSid.v){
    col.CT  <- ct
    if(ct=="allCT"){ col.CT <- CT.v } 
    col.nme <- paste("s_Cp_", rep(Cp.lst[[RS]], each=length(col.CT)), 
                     "_ct_", rep(col.CT, times=length(Cp.lst[[RS]])), 
                     "_e", sep="")
    temp <- as.matrix(BIN.MX[,col.nme])
    temp[ !temp%in%as.numeric(Cs.lst[[RS]]) ] <- 0
    RS.lst[[RS]] <- dimnames(temp)[[1]][rowSums(temp)>0]
    rm(temp, col.nme, col.CT)
  }
  #-------------------Prepare matrix for UpSet plot
  mx <- ComplexHeatmap::list_to_matrix(RS.lst)
  rm(RS.lst);gc()
  UPSOBJ <- ComplexHeatmap::make_comb_mat(mx, mode="distinct",
                                          remove_empty_comb_set=FALSE,
                                          remove_complement_set=FALSE,
                                          universal_set=uni.set)
  if(ct=="allCT"){
    save(x=UPSOBJ, file=paste0(out.dir, "/", out.name, ".RData"))
  }
  rm(mx); gc()
  
  # Collect UPSMX data
  CS[[ct]] <- comb_size(UPSOBJ)/totreg
  SS[[ct]] <- (set_size(UPSOBJ)[RSid.v])/totreg
  
  UPSOBJ <- UPSOBJ[comb_size(UPSOBJ)>0]
  cs <- comb_size(UPSOBJ)
  ss <- set_size(UPSOBJ)
  
  ht <- UpSet(UPSOBJ, 
              set_order=match(x=set_name(UPSOBJ), table=RSid.v), 
              comb_order=order(-comb_degree(UPSOBJ)),
              row_title=paste0(out.name, "\ndistinctmode\ntotalunionregion=", totreg),
              top_annotation=upset_top_annotation(UPSOBJ, ylim=ylim.v, height=unit(2, "cm"),
                                                  axis_param=list(side="left", at=at.v)
                                                  ),
              right_annotation=upset_right_annotation(UPSOBJ, ylim=ylim.v, width=unit(2, "cm"),
                                                      axis_param=list(side="bottom", at=at.v)
                                                      ),
              left_annotation=NULL, show_row_names=TRUE)
  ht <- draw(ht)
  c.od <- column_order(ht)
  r.od <- row_order(ht)
  decorate_annotation("Intersection\nsize", {
    grid.text(label=format(cs[c.od]/totreg, digits=3), x=1:length(cs), 
              y=unit(cs[c.od], "native") + unit(2, "pt"), 
              default.units="native", just=c("center", "bottom"), 
              gp=gpar(fontsize=5, col="black"), rot=0)
  })
  decorate_annotation("Set size", {
    grid.text(label=format(ss[r.od]/totreg, digits=3), y=length(ss):1, 
              x=unit(ss[r.od], "native") + unit(5, "pt"), 
              default.units="native", just=c("center", "bottom"), 
              gp=gpar(fontsize=5, col="black"), rot=90)
  })
  
  rm(UPSOBJ, cs, ss, c.od, r.od, ht, out.name); gc()
  print(paste0(ct, " done!"), quote=FALSE)
  
} # CT.v for loop end
dev.off()

UPSMX <- list(CS=do.call("rbind", CS), SS=do.call("rbind", SS),
              totunionreg=totreg)
save(UPSMX, file=paste0(out.dir, "/", gcb, "_allCT_RS_upsmx.RData"))

# rm(list=ls()); gc()



