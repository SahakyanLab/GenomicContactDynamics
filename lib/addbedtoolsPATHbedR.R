addbedtoolsPATHbedR <- function(PATH = "/Users/ltamon/prog/bedtools2/bin",
                                returnPATHorig = TRUE){
  
  library(bedr)
  
  # Add path containing bedtools to PATH (temporary)
  PATHorig <- Sys.getenv("PATH")
  if( !grepl(pattern=bedtoolsPATH, x=Sys.getenv("PATH")) ){
    Sys.setenv(PATH=paste(PATHorig, bedtoolsPATH, sep=":"))
  }
  # Confirm that bedtools can be accessed
  if(!check.binary("bedtools", verbose=TRUE)){
    stop("bedtools not found.")
  }
  
  if(returnPATHorig==TRUE){
    return(gsub(x=PATHorig, pattern=paste0(":", bedtoolsPATH), 
                replacement=""))
  }
  
}