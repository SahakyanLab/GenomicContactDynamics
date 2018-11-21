## FUNCTION ####################################################################
# This function takes the path and structuring of the genome/chromosome file   #
# (fasta) directory. It reads the genome/chromosome of a given, chr.id,        #
# identifier by creating a globally accessible fasta object. The object name   #
# consists of the supplied genome prefix, chr.id and fastafile.ending. In case #
# that object already exists in any of the current or parrent environments, the#
# function will do nothing. In both cases, the function will return a string   #
# message with return().                                                       #
# If remove.other.loads = TRUE (default is FALSE), then all the instances of   #
# objects in .GlobalEnv that start with genome.prefix will be deleted in case  #
# the specific chr.id object is not found, hence a new fasta file has to be    #
# loadec. This is useful while reading another chromosome, while removing the  #
# others via the same function.                                                #
# REQUIRES:                                                                    #
# LIB.TRANTOR = "/Users/alex/GIT/GITrepo/TrantoR"                              #
# source(paste0(LIB.TRANTOR, "/GEN_readfasta.R"))                              #
# source(paste0(LIB.TRANTOR, "/UTIL_readLinesFast.R"))                         #
# The underlying functions may require fread() from <data.table> library.      #
# USAGE EXAMPLE:                                                               #
# loadGenome(PATH.genome = "/Volumes/Data/Database/human_genome_unmasked_37.73",
#            genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
#            fastafile.ending = ".fa",
#            chr.id = 1, silent = FALSE,
#            remove.other.loads = FALSE)
################################################################################
loadGenome <- function(PATH.genome = "/Volumes/Data/Database/human_genome_unmasked_37.73",
                       genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
                       fastafile.ending = ".fa",
                       chr.id = 1, # any of c(1:22,"X","Y","MT")
                       silent = FALSE,
                       remove.other.loads = FALSE,
                       split = TRUE # sequence splitting while loading fasta
                       # If loaded genome exists (hence will not be freshly
                       # created by the current loadGenome call, <split> should
                       # macth the loaded genome split pattern.
                      ){

  genome.filename <- paste0(genome.prefix,chr.id,fastafile.ending)
  
  # exists(..., inherits=TRUE) looks in all parent environments as well.
  if(!exists(genome.filename, inherits=TRUE)){

    #----------------------
    if(remove.other.loads){
      # pos=1 is the first ebtry of search(), hence ".GlobalEnv".
      to.remove <- grep(genome.prefix, ls(pos=1L), value=TRUE)
      if(length(to.remove)!=0){
        rm(list=to.remove, pos=1L)
      }
    }
    #----------------------
    
    eval(parse(text=
               paste0(genome.filename,
                      " <<- readfasta('",
                      PATH.genome,"/",
                      genome.filename,
                      "', fastread=TRUE, split=split)")
    ))
    msg <- paste0("loadGenome: ",genome.filename," is now loaded.")

  } else {
    msg <- paste0("loadGenome: ",genome.filename," has already been loaded.")
  }
  
  if(silent){ return(msg) } else { print(msg, quote=FALSE) }
  
}
################################################################################
