################################################################################
# Combine key numbers per MUTBIN.DF
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
numbers.dir = paste0(wk.dir, "/out_mutCalcPerBin")
out.dir = paste0(wk.dir, "/out_checkNumbers")
### OTHER SETTINGS #############################################################
data.id = "donor_centric_PCAWG" # "CosmicNCV", "donor_centric_PCAWG"
src.id = "Hg19" # "Hg19" | "hg38ToHg19"
out.id = "comprehensivefinal"
sigEpLim.id.v = c("-1_0", "1000rawInf")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
file.v <- list.files(path=numbers.dir, pattern=data.id, full.names=F)
file.v <- file.v[grepl(x=file.v, pattern=".csv")]

x <- list()
for(file in file.v){
  
  x[[file]] <- read.csv(file=paste0(numbers.dir, "/", file), 
                          stringsAsFactors=F)
  
} # mut.id.v
x <- do.call("rbind.data.frame", c(x, stringsAsFactors=F))
x <- x[order(x$loc.id, x$sigEpLim.id, x$mut.id, x$SIG.id),]

#ggplot(data=x[x$sigEpLim.id!="nosampfilter",], aes(x=Nsamp)) +
ggplot(data=x[x$sigEpLim.id%in%sigEpLim.id.v,], aes(x=Nsamp)) +
  geom_density() +
  scale_x_continuous(breaks=seq(from=0, to=1800, by=25)) +
  labs(title=paste0(data.id, "_", src.id, "_", out.id, "_", 
                    'x[x$sigEpLim.id!="nosampfilter","Nsamp"]')) + 
  theme(axis.text.x=element_text(size=5, angle=90))

ggsave(filename=paste0(out.dir, "/", data.id, "_", src.id, "_", out.id, 
                       "_Nsamp_distribution.pdf"), height=10, width=10, unit="in")

write.csv(x=x, file=paste0(out.dir, "/", data.id, "_", src.id, "_", out.id, ".csv"),
          row.names=F)

# rm(list=ls()); gc()