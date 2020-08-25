################################################################################
# Use UpSetR plot to visualise overlap each region set across tissues. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
binmx.dir = paste0(wk.dir, "/binmx/out_bindata_1perc_HiCNorm")
out.dir = paste0(wk.dir, "/out_RSolapAcrossCT/CpAllCs1perc_HiCNorm")
### OTHER SETTINGS #############################################################
gcb = "min2Mb" 
chr.v = paste0("chr", c(1:22, "X"))
CT.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
         "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
Cp.v = 1:21
Cs.v = 0.01
id = "CpAllCs1perc" #"CpAll" #"CptopCP3" #"Cp1" #"CpAllCs1perc" #"CpAllCsmatch"

# If combsizethreshfr = 0.1, don't display top bar for combination comprised of equal 
# or less than 10% of total regions. Do this only if necessary for clarity. 
combsizethreshfr = 0

plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ComplexHeatmap)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.name <- paste0(gcb, "_", id)
print(paste0(out.name, "..."), quote=FALSE)

Cp.v <- sort(unique(Cp.v), decreasing=FALSE)
Cp.v.len <- length(Cp.v)
CT.v <- sort(unique(CT.v))
CT.v.len <- length(CT.v)

if(plotOnly==FALSE){
  
  # Prepare contacting region coordinates filtering by Cp.v
  if( is.null(chr.v) ){ chr.v <- paste0("chr", c(1:22, "X")) }
  BIN.MX <- sapply(X=chr.v, simplify=FALSE, FUN=function(chr){
    load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
    BIN.MX <- cbind.data.frame(chr=chr, BIN.MX, stringsAsFactors=FALSE)
  })
  BIN.MX <- do.call("rbind", BIN.MX)
  if( any(duplicated(dimnames(BIN.MX))) ){ stop("Duplicated region names.") }
  #---------------------------------------
  temp <- rep(CT.v, each=Cp.v.len)
  col.nme <- paste("s_Cp_", rep(Cp.v, times=CT.v.len), "_ct_", temp, 
                   "_e", sep="")
  #---------------------------------------
  # Prepare matrix for UpSet plot
  BIN.MX <- data.matrix(BIN.MX[,col.nme]); rm(col.nme)
  BIN.MX[!BIN.MX%in%as.numeric(Cs.v)] <- 0
  temp1 <- sapply(X=colnames(BIN.MX), USE.NAMES=FALSE, FUN=function(x){
    strsplit(x=x, split="_", fixed=TRUE)[[1]][5]
  })
  if( !identical(temp, temp1) ){ stop("Wrong order.") }
  mx <- sapply(X=CT.v, simplify=FALSE, FUN=function(ct){
    rsPerCT <- rowSums( x=as.matrix(BIN.MX[,temp==ct]) )
    rsPerCT[rsPerCT>0] <- 1
    return(rsPerCT)
  })
  mx <- do.call("cbind", mx)
  if( nrow(mx)!=nrow(BIN.MX) ){ stop("Number of regions different.") }
  mx <- mx[rowSums(mx)>0,] # Remove regions not in any tissue
  rm(temp, temp1, BIN.MX)
  
  upsobj <- ComplexHeatmap::make_comb_mat(mx, mode="distinct",
                                          remove_empty_comb_set=TRUE,
                                          # Removed by previous line
                                          remove_complement_set=FALSE)
  UPSOBJ <- list( upsobj=upsobj, totunionreg=nrow(mx) )
  rm(mx); gc()
  save(x=UPSOBJ, file=paste0(out.dir, "/", out.name, "_upsetplot.RData"))
  
} else {
  load(file=paste0(out.dir, "/", out.name, "_upsetplot.RData"))
}
#---------------------------------------UpSet plot
pdf(file=paste0(out.dir, "/", out.name, "_csthresh", 
                format(combsizethreshfr, scientific=FALSE), 
                "fr_upsetplot.pdf"), width=10, height=10)
totreg <- UPSOBJ$totunionreg
UPSOBJ <- UPSOBJ$upsobj

UPSOBJ <- UPSOBJ[(comb_size(UPSOBJ)/totreg)>combsizethreshfr]
ss <- set_size(UPSOBJ)
ht <- UpSet(UPSOBJ, set_order=order(-ss), comb_order=order(-comb_degree(UPSOBJ)),
            row_title=paste0(out.name, "\ndistinctmode\ntotalunionregion=",
                             totreg, "\nfiltercombsize>", combsizethreshfr, "fr"),
            left_annotation=NULL, show_row_names=TRUE)
ht <- draw(ht)
c.od <- column_order(ht)
cs <- comb_size(UPSOBJ)
decorate_annotation("Intersection\nsize", {
  grid.text(label=format(cs[c.od]/totreg, digits=3), x=seq_along(cs), 
            y=unit(cs[c.od], "native") + unit(2, "pt"), 
            default.units="native", just=c("center", "bottom"), 
            gp=gpar(fontsize=5, col="black"), rot=0)
})
dev.off()

# rm(list=ls()); gc()




