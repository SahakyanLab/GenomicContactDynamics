################################################################################
# Per chr, per norm type, generate plot of sd vs. mean/median Hi-Cs value of
# contact. Per stat (mean/median), 4 plots are generated based on the contacts 
# included.
# If(num.0s<=0), only contacts present across tissues are included. If(num.0s<=21),
# all contacts present in at least one tissue are included. There should be no
# contacts considered not present in any tissue (since we took the upper.tri of  
# MELT.MX) so the equal to 21 part is unnecessary. Cutoffs for num.0s tested are 
# 0, 5, 10, and 21. The code also saves an object R.MX, a list of matrices per 
# norm type containing Pearson's R (sd~mean/median) for each plot generated. 
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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/1_HiC_Human21_Expl"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/1_HiC_Human21_Expl"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
ijstat.dir = paste0(wk.dir, "/out_rawVsnorm")
out.dir = paste0(wk.dir, "/out_sdDivmean_plot")
### OTHER SETTINGS #############################################################
chr.v = c(paste0("chr", c(1:22, "X")))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(ggpubr)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
norm.v <- c("RAW_primary_cohort", "HiCNorm_primary_cohort", "HiCNorm_QQ_primary_cohort")
chr.v.len <- length(chr.v)
norm.v.len <- length(norm.v)

# Initialise matrix for R values
R.MX <- vector(mode="list", length=3)
R.MX[[1]] <- matrix(data=NA, nrow=length(chr.v), ncol=8, 
                    dimnames=list(chr.v, 
                                  paste0(rep(c("MEAN","MEDIAN"), times=4), rep(c(0,5,10,21), each=2)))
)
R.MX[[2]] <- R.MX[[3]] <- R.MX[[1]]
names(R.MX) <- norm.v

for(chr in chr.v){
  
  for(nrm in norm.v){
    
    flenme <- paste0(ijstat.dir, "/human_", chr, "_uppertri_", nrm, "_ijstat.RData")
    if( file.exists(flenme) ){
      load(file=paste0(ijstat.dir, "/human_", chr, "_uppertri_", nrm, "_ijstat.RData"))
    } else {
      print(paste0(flenme, " DNE."), quote=FALSE)
      next
    }; rm(flenme)
    ij.len <- nrow(IJSTAT.MX)
    row.num0s <- IJSTAT.MX$row.num0s
    IJSTAT.MX <- data.frame(sd=c(IJSTAT.MX$row.sds, IJSTAT.MX$row.sds),
                            value=c(IJSTAT.MX$row.means, IJSTAT.MX$row.medians),
                            stat=c(rep("MEAN", times=ij.len), rep("MEDIAN", times=ij.len)),
                            stringsAsFactors=FALSE); gc()
    
    max.y <- ceiling(max(IJSTAT.MX$sd))
    max.x <- ceiling(max(IJSTAT.MX$value))
    
    p.lst <- list()
    for( i in c(0,5,10,21) ){
      num0.TF <- row.num0s<=i
      i <- as.character(i)
      p.lst[[i]] <- ggplot(data=IJSTAT.MX[c(num0.TF,num0.TF),], aes(x=value, y=sd)) +
        geom_point(alpha=0.5, shape=1) + 
        scale_x_continuous(limits=c(0, max.x)) +
        scale_y_continuous(limits=c(0, max.y)) +
        labs(title=paste0("human_", chr, "_uppertrichr_", nrm, "_num.0s<=", i), 
             y="sd") +
        bgr2 +
        facet_grid(.~stat)
      
      # Get Pearson's correlation
      stat.ind <- (1:ij.len)[num0.TF]
      R.MX[[nrm]][chr,paste0("MEAN",i)] <- cor(x=IJSTAT.MX$value[stat.ind],
                                               y=IJSTAT.MX$sd[stat.ind])
      stat.ind <- stat.ind+ij.len
      R.MX[[nrm]][chr,paste0("MEDIAN",i)] <- cor(x=IJSTAT.MX$value[stat.ind],
                                                 y=IJSTAT.MX$sd[stat.ind])
      rm(num0.TF, stat.ind, i); gc()
    }
    rm(IJSTAT.MX, max.x, max.y); gc()
    p.arr <- ggarrange(plotlist=p.lst, nrow=4, ncol=1)
    ggexport(p.arr, width=8000, height=16000, res=800,
             filename=paste0(out.dir, "/human_", chr, "_uppertri_", nrm, "_ijstat.png"))
    rm(p.arr, p.lst)
    
  } # norm.v for loop end
  
  print(paste0(chr, " done!"), quote=FALSE)
  
} # chr.v for loop end

save(R.MX, file=paste0(out.dir, "/human_uppertri_RMX_ijstat.RData"))

# rm(list=ls()); gc()
