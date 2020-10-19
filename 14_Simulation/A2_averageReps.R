################################################################################
# Get average frequency of replicates per simulation type. Put replicates of 
# a simulation type in one directory then loop over directories to get average
# simulation matrix per type. Name of directory is used to name average matrix.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else if(whorunsit == "LiezelCluster"){
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_Simulation"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
simmap.dir = paste0(wk.dir, "/simulation_contact_maps")
out.dir = paste0(wk.dir, "/out_averageReps")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
dir.v <- list.files(simmap.dir)

for(dr in dir.v){
  
  print(paste0(dr, "..."), quote=FALSE)
  
  rep.v <- list.files(paste0(simmap.dir, "/", dr))
  
  REP <- sapply(X=rep.v, simplify=FALSE, FUN=function(rep){
    sim.mx <- data.matrix(fread(file=paste0(simmap.dir, "/", dr, "/", rep), 
                                header=FALSE, data.table=FALSE, stringsAsFactors=FALSE))
    dimnames(sim.mx) <- NULL
    return(sim.mx)
  })
  REP.len <- length(REP)
  print(paste0(REP.len, " replicates."), quote=FALSE)
  
  # Check if dimension of replicates the same
  dimval <- lapply(X=REP, FUN=dim)
  dimval <- unique(unlist(dimval))
  if( length(dimval)!=1 ){ stop("Problem with dimension of matrices.") }
  
  # Calculate average matrix of replicates 
  a <- paste(paste0("REP[[", 1:REP.len, "]]"),collapse="+")
  eval(parse(text=paste0(
    "ave.mx <-(", a, ")/", REP.len
  )))
  rm(REP, REP.len, a, rep.v, dimval); gc()
  
  write.table(ave.mx, file=paste0(out.dir, "/", dr, "_avemx.txt"), col.names=FALSE, 
              row.names=FALSE, quote=FALSE, sep="\t")
  rm(ave.mx); gc()
  
}

# rm(list=ls()); gc()
