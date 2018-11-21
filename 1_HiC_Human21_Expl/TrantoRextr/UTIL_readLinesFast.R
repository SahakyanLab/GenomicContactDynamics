################################################################################
## Faster version of readLines. Is dangerous to use in case the file to be    ##
## read is modified during this code execution.                               ##
## The first part is adapted from:                                            ##
## http://mlt-thinks.blogspot.co.uk/2011/08/faster-files-in-r.html            ##
################################################################################

readLinesFast <- function(filename, freadchar="@") {

  sz <- file.info( filename )$size

  if(.Machine$integer.max >= sz){

    return(
      strsplit(
        readChar( filename, nchars=sz, useBytes=TRUE ),
        "\n", fixed=TRUE, useBytes=TRUE
       )[[1]]
     )

  } else {

    # For large files, readChar has integer limitation problem. Using fread of
    # data.table library instead. Note: readr and tweaking readChar did not work.
    print("readLinesFast: the file is pretty big,",quote=F)
    print("      parsing via fread of data.table.", quote=F)
    print("     If the file contains the selected", quote=F)
    print(paste0("              <freadchar>=",freadchar," character,"),quote=F)
    print("there may be problems with fread parsing!", quote=F)
    # http://stackoverflow.com/questions/32920031/how-to-use-fread-as-readlines-without-auto-column-detection
    #library(data.table)
    return(data.table::fread(filename, sep=freadchar, header=F)[[1L]])

  }

}
################################################################################
