################################################################################
baseageFilepath = "BaseAgeEns76_Dataset/base_age_76.bed.gz"
splitdataPath   = "BaseAgeEns76_splitData"

library(data.table)
################################################################################

system( paste0("gunzip -k ./",baseageFilepath) )
ungzipped.filepath <- strsplit(baseageFilepath, split=".gz")[[1]]
data               <- fread(ungzipped.filepath, sep="\t", header=F)
file.remove(ungzipped.filepath)

dir.create(splitdataPath)

for(  chr in paste0("chr",c(1:22,"X","Y","MT"))  ){
  print(chr, quote=FALSE)
  write.table(data[data$V1==chr,],
              quote=TRUE, row.names=FALSE, col.names=FALSE,
              file=paste0(splitdataPath, "/",chr,"_BaseAgeEns76.txt"))
}

rm(data)
################################################################################
# Saving the unique Base Age categories from chr1 data alone:

chr1.data <- read.table(paste0(splitdataPath, "/chr1_BaseAgeEns76.txt"),
                        sep=" ", header=F, as.is=T)
unique.baseageinfo <- unique(
                      paste(chr1.data[,4], chr1.data[,5],chr1.data[,6], sep="|")
                      )

unique.baseageinfo <- read.table(textConnection(unique.baseageinfo), sep="|")
# Sort by baseAge score ascending order:
unique.baseageinfo <- unique.baseageinfo[order(unique.baseageinfo[,2]), c(2,3,1)]


# BaseAgeScore ColourRGB
write("BaseAgeScore ColourRGB Species", file="UniqueBaseAgeChr1.txt")
write.table(unique.baseageinfo,
            quote=TRUE, row.names=FALSE, col.names=FALSE,
            file="UniqueBaseAgeChr1.txt", append=TRUE)


