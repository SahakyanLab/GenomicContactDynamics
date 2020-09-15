################################################################################
# Inspect simulation maps from Zahra
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
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
simmap.dir = paste0(wk.dir, "/simulation_contact_maps")
out.dir = paste0(wk.dir, "/out_inspect_simMap")
### OTHER SETTINGS #############################################################
# Length should be equal to number of maps
col.v <- c("blueviolet", "darkred", "darkblue", "darkgreen", "darkgoldenrod", 
           "darkcyan")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
rep.v <- list.files(simmap.dir, full.names=TRUE, recursive=TRUE)
rep.v.len <- length(rep.v)
col.v <- adjustcolor(col.v, alpha=0.5)

pdf(file=paste0(out.dir, "/simvalues_dens.pdf"), width=10, height=10)

for(r in 1:rep.v.len){
  
  rep <- rep.v[r] 
  print(paste0(rep, "..."), quote=FALSE)
  
  sim.mx <- data.matrix(fread(file=rep, header=FALSE, data.table=FALSE, 
                              stringsAsFactors=FALSE))
  dimnames(sim.mx) <- NULL
  print(dim(sim.mx), quote=FALSE)
  
  lw.v <- sim.mx[lower.tri(sim.mx, diag=FALSE)]
  temp <- t(sim.mx)
  up.v <- temp[lower.tri(temp, diag=FALSE)]
  rm(temp); gc()
  
  stat.v <- c( max.nodiag=max(c(lw.v, up.v)), min.nodiag=min(c(lw.v, up.v)), 
               diag=unique(diag(sim.mx)) )
  print(stat.v, quote=FALSE)
  
  if(r==1){
    plot(density(log10(sim.mx[upper.tri(sim.mx, diag=FALSE)])), col=col.v[r],
         ylim=c(0,5))
    abline(v=log10(1.5))
  } else {
    lines(density(log10(sim.mx[upper.tri(sim.mx, diag=FALSE)])), col=col.v[r])
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
  
  rm(lw.v, up.v, rep, sim.mx); gc()
  
}

dev.off()

# rm(list=ls()); gc()
