################################################################################
#Title
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start_time <- Sys.time()

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelCluster"){
    # Main directory for the task
    objective.dir = "/t1-data/user/ltamon/ProjectBoard/_GenomicContactDynamics/3_AnnotationVsPersist"
    # Annotation file
    annoFile.dir   = "/t1-home/icbml/ltamon/Database/ucsc_tables"
    # HiC_Human21 persistent contact files
    persist.dir   = "/t1-home/icbml/ltamon/Database/HiC_features_GSE87112_RAWpc"
    os = "Linux"
  } else if(whorunsit == "AlexMac"){
    objective.dir = "./"
    annoFile.dir   = "/Volumes/Data/Database/ucsc_tables"
    persist.dir   = "/Volumes/Data/Database/HiC_features_GSE87112_RAWpc"
    os = "Mac"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
### OTHER SETTINGS #############################################################

################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################

################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################

end_time <- Sys.time()
end_time-start_time 