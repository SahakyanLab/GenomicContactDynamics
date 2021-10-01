################################################################################
# Define map directory (metric.dir) depending on metric
### FUNCTION ###################################################################
getmapdir <- function(metric='type of map; for non-simulation maps metric name 
                      should match source directory name; e.g. for metric=Cs.norm
                      Cs.norm.dir should be defined', 
                      simmap.dir='directory of simulation maps'){
  
  if( grepl(x=metric, pattern="SIM.", fixed=T) ){
    
    metric.dir <- paste0(simmap.dir, "/", metric)
  
  } else if( metric%in%c("Cs.raw", "Cs.norm", "Cp") | 
             grepl(x=metric, pattern="CII.cont.|CII.disc.") ){
    
    eval(parse(text=paste0(
      'metric.dir <- ', metric, '.dir'
    )))
    
  } else {
    stop("getmapdir(): Invalid metric argument.")
  }
  
  return(metric.dir)
  
}
################################################################################
