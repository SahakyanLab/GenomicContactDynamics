################################################################################
# Calculate per cluster, fraction of transposon sites (relative to total
# transposon site only) with age 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/17_RepeatAge")
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
repeatmx.dir = paste0(wk.dir, "/z_ignore_git/out_addToSummary")
out.dir = paste0(wk.dir, "/out_frTransposonWithAge")
### OTHER SETTINGS #############################################################
transposon.class = c("LINE", "SINE", "DNA", "LTR", "RC",
                     "LINE?", "SINE?", "DNA?", "LTR?",
                     "Other", "Unknown", "Unknown?")
not.transposon.class = c("Low_complexity", "RNA", "rRNA", "Satellite",
                         "scRNA", "Simple_repeat", "snRNA", "srpRNA", "tRNA")


rnk = "GiorPubl372rank" #c("Giordano364rank", "GiorPubl372rank", "Publrank")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# REPEAT.MX (base)
load(file=paste0(repeatmx.dir, "/hg19repeats_repName.RData"))
rownames(REPEAT.MX) <- NULL
REPEAT.MX <- REPEAT.MX[,c("cluster", "copyNumber", "repName", "repClass", "repFamily",
                          rnk)]

transposon.class.tmp <- setdiff(unique(REPEAT.MX$repClass), not.transposon.class)
if( !identical(sort(transposon.class.tmp), sort(transposon.class)) ){
  
  stop("Checkpoint 1.")
  rm(REPEAT.MX)
  
} 

is.rank <- REPEAT.MX[[rnk]] > 0
is.transposon <- REPEAT.MX$repClass %in% transposon.class

# Check if all subfam in ranking are transposons


REPEAT.MX$copyNumber <- as.numeric(as.character(REPEAT.MX$copyNumber))

clust.v <- sort(unique(REPEAT.MX$cluster))

FRTR <- list()
pdf(file=paste0(out.dir, "/hg19repeats_frRepeatSubfam_perCluster.pdf"), 
    height=10, width=10)
par(mfrow=c(2,1))

for(cl in clust.v){
  
  is.cl <- REPEAT.MX$cluster %in% cl
  
  # Proportion of sites per subfamily, per cluster
  
  cl.total <- sum(REPEAT.MX$copyNumber[is.cl])
  perc.sf <- REPEAT.MX$copyNumber[is.cl] / cl.total * 100
  cut.off <- unname(sort(perc.sf, decreasing=T)[10])
  labs <- paste0(REPEAT.MX$repName[is.cl], " - ", format(round(perc.sf, 4), scientific=F), "%")
  labs[as.numeric(perc.sf) < cut.off] <- NA
  
  pie(x=REPEAT.MX$copyNumber[is.cl], labels=labs, srt=-10,
      main=paste0(cl, "_frRepeatSubfam_perCluster_relTo_totalSitesPerCluster=", cl.total)) 
  
  # Proportion of transposon sites with age per cluster
  
  cl.total.transposon <- sum(REPEAT.MX$copyNumber[is.transposon & is.cl])
  
  FRTR[[cl]] <- c(TrCountInRank=sum(REPEAT.MX$copyNumber[is.cl & is.rank]),
                  totTrSites=cl.total.transposon,
                  totAllSites=cl.total)
}

dev.off()

FRTR <- do.call("rbind", FRTR)

write.csv(FRTR, file=paste0(out.dir, "/hg19repeats_TransposonCount_perClusterIn", rnk, ".csv"),
          row.names=T)

# rm(list=ls()); gc()
