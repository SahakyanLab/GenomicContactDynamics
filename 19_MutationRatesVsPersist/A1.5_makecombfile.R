################################################################################
# Make a list of combinations of loc, mut, SIG and sigEpLim.
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
out.dir = paste0(wk.dir, "/out_makecombfile")
sigExposure.dir = paste0(data.dir, "/signal_mutSig/out_samplesForSignature")
### OTHER SETTINGS #############################################################
out.id = "-1_0" # "nosampfilter_1000rawInf"
loc.v = c("exon", "intron", "intergenic", "intron_intergenic", 
          "exon_intron_intergenic_intergenic.excl")

mut.v = c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

sigEpLim.v = "-1_0" #"1000rawInf" 

sigExpPERC.df <- read.csv(file=paste0(sigExposure.dir, "/donorlist_signatureExposurePERCENT.csv"), 
                          header=T, stringsAsFactors=F)
sigExpPERC.df <- sigExpPERC.df[,-c(1,3:7)]
print(paste0("Signature exposure dfs has ", ncol(sigExpPERC.df[,-1]), " signatures..."), quote=FALSE)
# Define signatures
SIG.v = c("RefSig.MMR1_RefSig.MMR2", 
          colnames(sigExpPERC.df)[colnames(sigExpPERC.df)!="alt.ID"])
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# nosampfilter
comb1 <- expand.grid(mut=mut.v, loc=loc.v, stringsAsFactors=F)
comb1 <- cbind(SIG="RefSig.1", comb1, sigEpLim="nosampfilter")
comb1 <- comb1[order(comb1[,"loc"], comb1[,"mut"]),]
# Rest
comb2 <- expand.grid(SIG=SIG.v, mut=mut.v, loc=loc.v, sigEpLim=sigEpLim.v,
                     stringsAsFactors=F)
comb2 <- comb2[order(comb2[,"loc"], comb2[,"sigEpLim"], comb2[,"mut"], comb2[,"SIG"]),]
# Combine comb1 and comb2
comb <- rbind.data.frame(comb1, comb2, stringsAsFactors=F)

comb.v <- paste(comb[,"SIG"], comb[,"mut"], comb[,"sigEpLim"], comb[,"loc"], sep=";")

write(x=comb.v, file=paste0(out.dir, "/combfile_", out.id))

# rm(list=ls()); gc()