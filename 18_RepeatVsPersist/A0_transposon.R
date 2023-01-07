################################################################################
# Make list (.RData) of families and subfamilies that are and are not transposons.
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
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")

repeatmx.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/17_RepeatAge/z_ignore_git/out_addToSummary")
out.dir = paste0(wk.dir, "/out_transposon")
### OTHER SETTINGS #############################################################
transposon.class = c("LINE", "SINE", "DNA", "LTR", "RC",
                     "LINE?", "SINE?", "DNA?", "LTR?",
                     "Other", "Unknown", "Unknown?")
not.transposon.class = c("Low_complexity", "RNA", "rRNA", "Satellite",
                         "scRNA", "Simple_repeat", "snRNA", "srpRNA", "tRNA")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(file=paste0(repeatmx.dir, "/hg19repeats_repName.RData"))
rownames(REPEAT.MX) <- NULL
REPEAT.MX <- REPEAT.MX[,c("repName", "repClass", "repFamily")]

transposon.class.tmp <- setdiff(unique(REPEAT.MX$repClass), not.transposon.class)
if( !identical(sort(transposon.class.tmp), sort(transposon.class)) ){
  
  stop("Checkpoint 1.")
  rm(REPEAT.MX)
  
} 

is.transposon <- REPEAT.MX$repClass %in% transposon.class

TRANSPOSON <- list()

TRANSPOSON[["yes"]] <- list(fam=REPEAT.MX$repFamily[is.transposon],
                            subfam=REPEAT.MX$repName[is.transposon],
                            subfamALL=REPEAT.MX$repName[is.transposon])
TRANSPOSON[["no"]] <- list(fam=REPEAT.MX$repFamily[!is.transposon],
                           subfam=REPEAT.MX$repName[!is.transposon],
                           subfamALL=REPEAT.MX$repName[!is.transposon])

save(TRANSPOSON, file=paste0(out.dir, "/hg19repeats_Transposon_yes_no.RData"))

# rm(list=ls()); gc()