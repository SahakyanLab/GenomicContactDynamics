################################################################################
# Determine final list of foi
################################################################################
################################################################################
finaliseFOI <- function(foi.dir = "/dir",
                        foifile = "path of foifile"){
  
  # List of features in directory
  foi.dir.v <- list.files(path=foi.dir, recursive=TRUE)
  foi.dir.v <- foi.dir.v[
    !grepl(x=foi.dir.v, pattern="\\.gz|ignore") & grepl(x=foi.dir.v, 
                                                        pattern="ct.*foi.*desc")
    ]
  
  if( !is.null(foifile) ){
    #foi.dir.v <- foi.dir.v[ grepl(x=foi.dir.v, 
    #                      pattern=paste(readLines(foifile), collapse="|")) ]
    foifile.v <- readLines(foifile)
    if( any(!foifile.v%in%foi.dir.v) ){
      stop("Feature in foifile not found in directory.")
    }
    foi.dir.v <- foifile.v 
  }
  
  return(foi.dir.v)
  
}