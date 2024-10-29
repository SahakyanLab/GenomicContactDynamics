################################################################################
# Summary plot showing contact-wise enrichment of all tested features
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar = TRUE) 
options(warn = 1) 

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = file.path(home.dir, "git/GenomicContactDynamics/lib/")
wk.dir = file.path(home.dir, "git/GenomicContactDynamics/7_FeaturePermutation")
featvscp.dir = paste0(wk.dir, "/z_ignore_git/out_foiVsij")
foifile.dir = paste0(wk.dir, "/foifile/foivsij")
out.dir = paste0(wk.dir, "/z_ignore_git/out_foivsij_organise")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chrs = paste0("chr", c("X", 1:22))
foifile_range = c(1:34, 100:200)
foivsij_plot_suffix = "any_querybin_foiVsij_ij"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
source(file.path(lib, "GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
foifile_lst <- sapply(as.character(foifile_range), USE.NAMES = TRUE, simplify = FALSE, function(i) {
  readLines(file.path(foifile.dir, paste0("foifile", i)))
})

foivsij_plot_names <- list.files(featvscp.dir, pattern = paste0(foivsij_plot_suffix, ".pdf"), recursive = FALSE)
if (any(duplicated(foivsij_plot_names))) {
  stop("Duplicated plots")
}

for (foi_char in names(foifile_lst)) {
  
  split_lst <- strsplit(foifile_lst[[foi_char]], split = "ct_|_foi_|_desc_")
  split_lst <- lapply(split_lst, function(x) x[x != ""])
  
  plot_name <- grep(paste0("_", foi_char, "_", foivsij_plot_suffix), foivsij_plot_names, value = TRUE)
  
  dir_name <- unique(unlist(lapply(split_lst, function(x) x[3])))
  
  if (length(dir_name) >= 1) {
    
    if (length(dir_name) > 1) {
      
      warning("foifile", foi_char, ": Multiple desc (i.e. dir_name). Combining desc as dir_name.")
      dir_name <- paste(dir_name, collapse = "-")
      
    }
      
    plot_out_dir <- file.path(out.dir, dir_name)
    if (!dir.exists(plot_out_dir)) {
      dir.create(plot_out_dir)
    }

  } else {
    stop("foifile", foi_char, ": No desc (i.e. dir_name)")
  }
  
  if (length(plot_name) == 1) {
    
    plot_name_new <- gsub(paste0("_", foi_char, "_", foivsij_plot_suffix),
                          paste0("_", dir_name, "_", foivsij_plot_suffix),
                          plot_name)
    
    file.copy(file.path(featvscp.dir, plot_name), file.path(plot_out_dir, plot_name_new), overwrite = TRUE)
    
  } else if (length(plot_name) == 0) {
    
    warning("foifile", foi_char, ": No existing plot found. Skipping.")
    next
    
  } else {
    stop("foifile", foi_char, ": Multiple plots returned.")
  }
  
  message("foifile", foi_char, " done!")
  
}

