## FUNCTION ####################################################################
# This function returns the required genomic sequence segment, fully exploiting#
# the capabilities of the loadGenome() function that skips unnecessarily loa-  #
# ding genomes twice.                                                          #
# REQUIRES:                                                                    #
# LIB.TRANTOR = "/Users/alex/GIT/GITrepo/TrantoR"                              #
# source(paste0(LIB.TRANTOR, "/GEN_loadGenome.R"))                             #
# source(paste0(LIB.TRANTOR, "/GEN_readfasta.R"))                              #
# source(paste0(LIB.TRANTOR, "/UTIL_readLinesFast.R"))                         #
# The underlying functions may require fread() from <data.table> library.      #
# USAGE EXAMPLE:                                                               #
#getGenomicSeq(PATH.genome = "/Volumes/Data/Database/human_genome_unmasked_37.73",
#              genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
#              fastafile.ending = ".fa",
#              chr.id = 18, # any of c(1:19,"X","Y","MT")
#              remove.other.loads = TRUE,
#              silent = FALSE, split = TRUE,
#              borders = c(1750000,1750500) )
################################################################################
getGenomicSeq <- function(PATH.genome = "/Volumes/Data/Database/human_genome_unmasked_37.73",
                          genome.prefix = "Homo_sapiens.GRCh37.73.dna.chromosome.",
                          fastafile.ending = ".fa",
                          chr.id = 1, # any of c(1:22,"X","Y","MT")
                          remove.other.loads = TRUE,
                          silent = FALSE,
                          split = TRUE, # sequence splitting while loading fasta
                          # If loaded genome exists (hence will not be freshly
                          # created by the current loadGenome call, <split> should
                          # macth the loaded genome split pattern.
                          borders = c(1750000,1750500)
                         ){

  genome.filename <- paste0(genome.prefix,chr.id,fastafile.ending)

  loadGenome(PATH.genome = PATH.genome,
             genome.prefix = genome.prefix,
             fastafile.ending = fastafile.ending,
             chr.id = chr.id, silent = silent,
             remove.other.loads = remove.other.loads, split=split)
  
  #---------
  if(split){
    return(eval(parse(text=
                    paste0(genome.filename,"$seq[borders[1]:borders[2]]")
          )))
  } else {
    return(eval(parse(text=
            paste0("substr(",genome.filename,
                   "$seq, start=borders[1], stop=borders[2])")
          )))    
  }
  #---------

}
################################################################################
