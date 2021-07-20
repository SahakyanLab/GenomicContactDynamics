################################################################################
# Get citation and version of R packages used
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/Rpkg"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
out.dir = paste0(wk.dir, "/out_Rlib_citationAndVersion")
### OTHER SETTINGS #############################################################
pkgPath = paste0(wk.dir, "/Rpkg.csv")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
out.id = "supp_gcd"
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
pkg.v <- sort(unique( as.character(read.csv(file=pkgPath, header=TRUE)[,1]) ), 
              decreasing=F) 
  
ip.mx <- utils::installed.packages()

# Packages not installed
nipkg.v <- pkg.v[!pkg.v%in%dimnames(ip.mx)[[1]]]
if( length(nipkg.v)>0 ){
  stop(paste0("Packages not installed: ", paste(nipkg.v, collapse=";")))
}
#install.packages(nipkg.v)

ver <- list()
for(p in pkg.v){
  
  print(paste0(p, "..."), quote=F)
  
  # Citation
  write(x=toBibtex(citation(p)), file=paste0(out.dir, "/", out.id, "_cite_Rpkg.bib"), 
        append=T)
  
  # Version
  ver[[p]] <- ip.mx[p,"Version"]
  
}

ver <- stack(ver)
ver$values <- as.character(ver$values)
ver$ind <- as.character(ver$ind)
ver <- ver[order(ver$ind, decreasing=F),]

# Output as latex string; remove backslash (\texttt), just add it on latex file
ver$ind <- paste0("TEXTTT{", ver$ind, "}")
ver$values <- paste0("(", ver$values, ")")

ver.str <- paste(ver$ind, ver$values, sep=" ")
ver.str <- paste(ver.str, collapse=", ")

write(ver.str, file=paste0(out.dir, "/", out.id, "_ver_tex_Rpkg"), 
      append=F)

# rm(list=ls()); gc()