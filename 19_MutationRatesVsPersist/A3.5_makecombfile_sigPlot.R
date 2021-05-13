################################################################################
# Make a list of combinations of gcb, data.id, src.id, calc, aggregatefunx, wise,
# sigEpLim.id, Nsamp.minfor the signature plots.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
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
out.dir = paste0(wk.dir, "/out_makecombfile_sigPlot")
### OTHER SETTINGS #############################################################
#gcb.v = "min2Mb"
#data.id.v = "donor_centric_PCAWG" # "CosmicNCV", "donor_centric_PCAWG"
#src.id.v = "Hg19" # "Hg19" | "hg38ToHg19"
#calc.v = "Nmsitenorm" # "Tmut" | "Nmsite" | "TmutDIVNmsite" | "Nmsitenorm" | "numWTSEQ"
#aggregatefunx.v = c("mean", "median")
wise.v = c("b", "cmean", "csd")
sigEpLim.id.v = c("-1_0", "1000rawInf")
Nsamp.min.v = 50
out.id = "min2Mb_donor_centric_PCAWG_Hg19_Nmsitenorm_mean"
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
comb <- expand.grid(#gcb=gcb.v, data.id=data.id.v, src.id=src.id.v, 
                    #calc=calc.v, aggregatefunx=aggregatefunx.v,
                    wise=wise.v, sigEpLim.id=sigEpLim.id.v, Nsamp.min=Nsamp.min.v,
                    stringsAsFactors=F)
comb <- comb[order(#comb$gcb, comb$data.id, comb$src.id, comb$calc, 
                   comb$Nsamp.min, comb$wise, 
                   #comb$aggregatefunx, 
                   comb$sigEpLim.id),]
comb <- paste(#comb$gcb, comb$data.id, comb$src.id, comb$calc, 
              comb$wise, 
              #comb$aggregatefunx,
              comb$sigEpLim.id, comb$Nsamp.min, sep=";")

write(x=comb, file=paste0(out.dir, "/combfile_", out.id))

# rm(list=ls()); gc()