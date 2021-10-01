################################################################################
# Inspect simulation map from Zahra
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
simmap.dir = paste0(wk.dir, "/Combined Matrices_Oxf/interacting_arms")
out.dir = paste0(wk.dir, "/out_inspect_simMap")
### OTHER SETTINGS #############################################################
out.id = "interacting_arms"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
library(RColorBrewer)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
fle.v <- list.files(path=simmap.dir, recursive=T)
fle.v.len <- length(fle.v)

# Colors of line plots, Length should be equal to number of maps
col.v <- RColorBrewer::brewer.pal(n=fle.v.len, name="Set1")

if( fle.v.len<length(col.v) ){
  stop("Colours should be equal to number of maps.")
}

pdf(file=paste0(out.dir, "/", out.id, "_simvalues_dens.pdf"), width=10, height=10)

for(f in 1:fle.v.len){
  
  print(paste0(out.id, ": ", fle.v[f], "..."), quote=F)
  
  sim.mx <- data.matrix(fread(file=paste0(simmap.dir, "/", fle.v[f]), header=F, 
                              data.table=F, stringsAsFactors=F))
  dimnames(sim.mx) <- NULL
  print(dim(sim.mx), quote=F)
  
  lw.v <- sim.mx[lower.tri(sim.mx, diag=FALSE)]
  temp <- t(sim.mx)
  up.v <- temp[lower.tri(temp, diag=FALSE)]
  rm(temp); gc()
  
  stat.v <- c( max.nodiag=max(c(lw.v, up.v)), min.nodiag=min(c(lw.v, up.v)), 
               diag=unique(diag(sim.mx)) )
  print(stat.v, quote=FALSE)
  
  if(f==1){
    
    plot(density(log10(sim.mx[upper.tri(sim.mx, diag=FALSE)])), col=col.v[f],
         ylim=c(0,20), xlim=c(0,5))
  
  } else {
    lines(density(log10(sim.mx[upper.tri(sim.mx, diag=FALSE)])), col=col.v[f])
  }

  if( !isSymmetric(sim.mx) ){
    
    ind.ns <- which((lw.v!=up.v))
    print(paste0(length(ind.ns), " contacts not symmetric."), quote=FALSE)
    
    mx.diff <- max(lw.v[ind.ns]-up.v[ind.ns])
    print(paste0("Max difference in value of asymmetric contacts is ", 
                 mx.diff), quote=FALSE)
    v <- c(lw.v[ind.ns], up.v[ind.ns])
    v <- c(max.val.asy=max(v), min.val.asy=min(v))
    print(v, quote=FALSE)
    
    rm(ind.ns, mx.diff, v); gc()
    
  } else {
    print("Symmetric.", quote=FALSE)
  }
  
  rm(lw.v, up.v, sim.mx); gc()
  
}

dev.off()

# rm(list=ls()); gc()
