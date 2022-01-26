################################################################################
# Calculate lengths from bed files
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=TRUE)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/16_GeneVsPersist")
    os = "Mac"
  #} else if(whorunsit == "LiezelCluster"){
  #  home.dir = "/project/sahakyanlab/ltamon" 
  # wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/3_AnnotationVsPersist")
  # os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
anno.nme = "hg19anno"
lib = paste0(home.dir, "/DPhil/lib")
data.dir = paste0(home.dir, "/Database")
out.dir = paste0(wk.dir, "/out_geneLength")
bed.dir = paste0(out.dir, "/beds_toGetRepFree_zerobased")
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
bednme.v = c(
  "ct_hg19_foi_Transcript_full_desc_DNA",
  "ct_hg19_foi_Transcript_full_repfree_desc_DNA",
  "ct_hg19_foi_Transcript_exon_desc_DNA",
  "ct_hg19_foi_Transcript_exon_repfree_desc_DNA",
  "ct_hg19_foi_Transcript_intron_desc_DNA",
  "ct_hg19_foi_Transcript_intron_repfree_desc_DNA"
)
bednme.v.len <- length(bednme.v)

bed.id.v <- strsplit(x=bednme.v,  split='_foi_Transcript_|_desc_')
bed.id.v <- unlist( lapply(X=bed.id.v, FUN=function(id) id[[2]] ) )

for(b in 1:bednme.v.len){
  
  bed.id <- bed.id.v[[b]]
  
  bed <- data.table::fread(file=paste0(bed.dir, "/", bednme.v[b]), header=F, 
                           data.table=F, stringsAsFactors=F)
  setnames(x=bed, old="V5", new="uniqueID")
  
  bed$len <- bed$V3-bed$V2 # 0-coordinates so no need for plus 1 to get length
  if( any(is.na(bed$len) | bed$len==0) ){
    stop(paste0(bed.id, ": Invalid length."))
  }
  
  # Total length and number of regions per unique id
  len.df <- aggregate(x=bed$len, by=list(bed$uniqueID), simplify=T, FUN=function(regions){
    c(sum(regions, na.rm=T), 
      length(regions))
  })
  len.df <- data.frame(uniqueID=len.df$Group.1, len.df$x, stringsAsFactors=F)
  colnames(len.df) <- c("uniqueID", 
                         paste0("len_", bed.id), 
                         paste0("num_", bed.id)) 
  
  if(b==1){
    LEN.DF <- len.df
  } else {
    
    if( !identical(sort(LEN.DF$uniqueID), sort(len.df$uniqueID)) ){
      warning(paste0(bed.id, ": Unique ids different than LEN.DF."))
    }
    LEN.DF <- merge(x=LEN.DF, y=len.df, by="uniqueID", all.x=T, all.y=T)
    
  }
  
}

# This is expected because I verified that transcript start and end are equal 
# to start of most upstream exon and end of most downstream exon, respectively. 
if( any(LEN.DF$num_exon != (LEN.DF$num_intron  + 1L), na.rm=T) ){
  stop("Intron number not equal to exon number minus 1.")
}

# Calculate absolute and relative repeat content
for( part in c("full", "exon", "intron") ){
  
  LEN.DF[[ paste0("len_repeat_", part) ]] <-
    LEN.DF[[ paste0("len_", part) ]] - LEN.DF[[ paste0("len_", part, "_repfree") ]]
  
  if( any(LEN.DF[[ paste0("len_repeat_", part) ]] < 0, na.rm=T) ){
    stop(paste0(part, ": repfree > orig length."))
  }
  
  LEN.DF[[ paste0("fr_repeat", "_", part) ]] <- 
    LEN.DF[[ paste0("len_repeat_", part) ]] / LEN.DF[[ paste0("len_", part) ]]
    
}

LEN.DF <- LEN.DF[,sort(colnames(LEN.DF))]
save(LEN.DF, file=paste0(out.dir, "/", anno.nme, "_ALL_annoLengths_base.RData"))

#Warning messages:
#1: exon_repfree: Unique ids different than LEN.DF. 
#2: intron: Unique ids different than LEN.DF. 
#3: intron_repfree: Unique ids different than LEN.DF. 

# rm(list=ls()); gc()