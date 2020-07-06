################################################################################
# Plot of permutation test result
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
### OTHER SETTINGS #############################################################
wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/out_associate_runA_TOP2B"
out.dir = paste0(wk.dir, "/plot_permtest")
out.name = "min2Mb_hg19_TOP2B_TF_MCF7_nperm10000_cp21_seed342"
evalf = "total length of intersection"
evalf.lab = "comOlap"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
### FUNCTION ###################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(paste0(wk.dir, "/", out.name, "_permtest.RData"))

pdf(file=paste0(out.dir, "/", out.name, "_", evalf.lab, ".pdf"), width=10, height=10)
plot(PERMT[[evalf]])
dev.off()

# rm(list=ls()); gc()