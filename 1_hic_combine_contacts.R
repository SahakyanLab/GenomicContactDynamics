################################################################################
# Path of the file holding the structure of the GSE87112 contact maps:
data.str.path       = "/Users/alex/GIT/GITrepo/HiC_Human21/data_structure.txt"
# The directory holding the contact maps of interest from GSE87112:
contact.map.dir     = "/Volumes/Data/Database/GSE87112/contact_maps/HiCNorm_QQ"
# Mid- and end-parts of the matrix files in the above directories:
midpart.mx.filename = ".nor.chr"
endpart.mx.filename = ".qq.mat"
# Path to the location, where the outputs are to be saved (make sure that the
# location exists):
out.path = "/Volumes/Data/Database/GSE87112/combined_contacts/HiCNorm_QQ_primary_cohort"
################################################################################
library(reshape2)
# Also directly calls fread() function from data.table library.
################################################################################




################################################################################

data.str <- read.table(data.str.path, header=TRUE, as.is=TRUE)
data.str <- data.str[data.str$source=="primary_cohort", ]

chrs <- c(1:22,"X")
for(chr in chrs){
  counter <- 1
  for(src in 1:21){
    label <- data.str[src,"label"]
    matrix.filepath <- paste0(contact.map.dir,"/primary_cohort/",
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
    #-----------------

    print(paste0("chr",chr," - ",label," is done."),quote=FALSE)
    rm(mx,melt.mx,na.rows)
    gc()
    counter <- counter + 1
  }

  contact.sums       <- rowSums(MELT.MX$upper.tri[,3:dim(MELT.MX$upper.tri)[2]])
  all.zero.rows      <- which(contact.sums==0)
  MELT.MX$upper.tri.nocontact <- MELT.MX$upper.tri[all.zero.rows,1:3]
  dimnames(MELT.MX$upper.tri.nocontact)[[2]] <- c("i","j","value")
  MELT.MX$upper.tri  <- MELT.MX$upper.tri[-all.zero.rows,]
  save(MELT.MX, file=paste0(out.path,"/human_chr",chr,"_allcontacts.RData"))
  rm(MELT.MX)
  gc()
}

################################################################################
