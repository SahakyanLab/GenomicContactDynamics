################################################################################
# Combine contact data from all tissues.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/1_HiC_Human21_Expl")
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/1_HiC_Human21_Expl")
  } else if(whorunsit == "LiezelLinuxDesk"){
    home.dir = "/home/ltamon"
    wk.dir = paste0(home.dir, "/DPhil/GCD_polished/1_HiC_Human21_Expl")
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(home.dir, "/Database")

bin.len.id = "50kb"
mx.type.id = "KRoe"

# Path of the file holding the structure of the GSE87112 contact maps:
data.str.path = paste0(wk.dir, "/data_structure.txt") 
## The directory holding the contact maps of interest from GSE87112:
#contact.map.dir = paste0(data.dir, "/GSE87112/contact_maps/HiCNorm_QQ")
contact.map.dir = paste0(data.dir, "/human_hg38_contacts/contact_maps_", bin.len.id, "/", mx.type.id)
# Mid- and end-parts of the matrix files in the above directories:
midpart.mx.filename = ".nor.chr"
endpart.mx.filename = paste0(".", mx.type.id, ".mat") #".qq.mat"
# Path to the location, where the outputs are to be saved (make sure that the
# location exists):
#out.path = paste0(data.dir, "/GSE87112/combined_contacts/HiCNorm_QQ_primary_cohort")
out.path = paste0(data.dir, "/human_hg38_contacts/combined_contacts_", bin.len.id, "/", mx.type.id)

source.id = "hg38_contacts" # "primary_cohort"
species.id = "human" #"human" "dme"
### OTHER SETTINGS #############################################################
chrs = c(1:22, "X") #as.character(c("X", "2L", "2R", "3L", "3R", "4"))
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
data.str <- read.table(data.str.path, header=TRUE, stringsAsFactors=FALSE)
data.str <- data.str[data.str$source==source.id, ]
src.len <- length(data.str[,1])

if( !dir.exists(out.path) ){
  stop("out.path does not exist.")
}

for(chr in chrs){
  counter <- 1
  for(src in 1:src.len){
    label <- data.str[src,"label"]
    matrix.filepath <- paste0(contact.map.dir,"/", #primary_cohort/",
                              data.str[src,"cell_type"],
                              midpart.mx.filename,chr,endpart.mx.filename)
    mx <- as.matrix(data.table::fread(matrix.filepath, sep="\t", header=F))
    dimnames(mx) <- list(1:dim(mx)[1], 1:dim(mx)[1])
    
    # Testing that the read matrix is symmetric.
    if(isSymmetric(mx)==FALSE){
      stop("The read matrix is not symmetric!")
    }
    
    # Setting the diagonal and lower triangular part to NA (hence for right
    # query, j > i should be kept).
    mx[!upper.tri(mx)] <- NA
    melt.mx <- melt(mx)
    dimnames(melt.mx)[[2]] <- c("i","j",label)
    na.rows <- which(is.na(melt.mx[,3]))
    
    #-----------------
    if(counter == 1) {
      
      MELT.MX <- NULL
      MELT.MX$upper.tri <- melt.mx[-na.rows,]
      MELT.MX$rest      <- melt.mx[na.rows,]
      dimnames(MELT.MX$rest)[[2]] <- c("i","j","value")
      
    } else {
      
      melt.mx <- melt.mx[-na.rows,]
      if(all(melt.mx$i==MELT.MX$upper.tri$i) &
         all(melt.mx$j==MELT.MX$upper.tri$j) ){
        eval(parse(text=paste0("MELT.MX$upper.tri <- cbind(MELT.MX$upper.tri, ",
                               label, "=melt.mx[,3])")))
      } else {
        stop("Matrices do not match!")
      }
      
    }
    if( length(MELT.MX$upper.tri[[label]])==0 ){ stop(label) }
    #-----------------
    
    print(paste0("chr",chr," - ",label," is done."),quote=FALSE)
    rm(mx,melt.mx,na.rows)
    gc()
    counter <- counter + 1
  }
  
  contact.sums       <- rowSums(MELT.MX$upper.tri[,3:dim(MELT.MX$upper.tri)[2], drop=FALSE])
  all.zero.rows      <- which(contact.sums==0)
  MELT.MX$upper.tri.nocontact <- MELT.MX$upper.tri[all.zero.rows,1:3]
  dimnames(MELT.MX$upper.tri.nocontact)[[2]] <- c("i","j","value")
  if(length(all.zero.rows)!=0){
    MELT.MX$upper.tri  <- MELT.MX$upper.tri[-all.zero.rows,]
  }
  save(MELT.MX, file=paste0(out.path,"/", species.id, "_chr", chr, "_allcontacts.RData"))
  rm(MELT.MX); gc()
  
}

# rm(list=ls()); gc()
