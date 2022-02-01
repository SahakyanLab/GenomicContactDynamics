################################################################################
# Anatogram for different organisms showing organs/parts of interest
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Set recommended global options

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=T)

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/misc/anatogram")
out.dir = paste0(wk.dir, "/out_make_anatogram")
### OTHER SETTINGS #############################################################
organs = readLines(con=paste0(wk.dir, "/organs"))
out.name = "female_organs" #"female_organs_perOrgan"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#devtools::install_github("jespermaag/gganatogram")
library(gganatogram)
library(ggsci)
library(yarrr)
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
data("hgFemale_key")

if( any(!organs%in%hgFemale_key$organ) ){
  stop("Some organ/s not in package.")
} else {
  organ.df <- hgFemale_key[hgFemale_key$organ %in% organs,]
}

rownames(organ.df) <- organ.df$organ
organ.df <- organ.df[organs,]
  
organ.df$value <- NULL

organ.df$type <- organ.df$organ
type.len <- length(unique(organ.df$type))

organ.df$colour <- colorRampPalette(yarrr::piratepal("pony", trans=0.5))(type.len) #colorRampPalette(brewer.pal(n=12, name="Set3"))(type.len)
organ.df$colour <- adjustcolor(col=organ.df$colour, alpha.f=0.5)

p <- gganatogram(data=organ.df, fillOutline="gray95", organism="human", 
                 sex="female", fill="colour") + 
  #facet_wrap(~type) + 
  theme_void()

ggsave(filename=paste0(out.dir, "/", out.name, ".pdf"), plot=p,
       height=10, width=5, units="in")


# rm(list=ls()); gc()