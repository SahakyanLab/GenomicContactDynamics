################################################################################
# Table of Tmut, Nmsite, TmutDIVNmsite and Nmsitenorm per bin. Indicate also which
# Cp/s the bin belongs to. Using this table make boxplots of the 4 calculations in C1 
# per Cp.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist1"
    binmx.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/binmx/out_bindata_1perc_HiCNorm"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
mutbin.dir = paste0(wk.dir, "/out_mutCalcPerBin")
out.dir = paste0(wk.dir, "/out_mutCalcPerCp")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
data.id = "donor_centric_PCAWG_sigEperc50_MMR12" # "CosmicNCV" | "donor_centric_PCAWG_as_is" | "donor_centric_PCAWG_sigEperc30_MMR12"
src.id = "Hg19" # "Hg19" | "hg38ToHg19"
Cp.v = 1:21
mut.v = c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
calc.v = c("Tmut", "Nmsite", "TmutDIVNmsite", "Nmsitenorm")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(reshape2)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
pdf(file=paste0(out.dir, "/", gcb, "_", data.id, "_", src.id, "_boxplots.pdf"),
    width=35, height=10)
par(mfcol=c(4,7))

for(mut in mut.v){
  
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=TRUE)
  # Load MUTBIN.DF
  load(file=paste0(mutbin.dir, "/", data.id, "_", src.id, "_", mut.id, "_mutCalcPerBin.RData"))
  
  MUTCP.DF <- list()
  chr.v <- unique(MUTBIN.DF$chr)
  for(chr in chr.v){
    
    load(file=paste0(binmx.dir, "/", chr, "_", gcb, "_bindata.RData"))
    #BIN.MX <- BIN.MX[rowSums(BIN.MX)>0,-(1:2)]
    BIN.MX <- BIN.MX[,-(1:2)]
    # c("totmut", "totmutdivNmsite", "mutbin"
    MUTCP.DF[[chr]] <- matrix(data=0, nrow=nrow(BIN.MX), ncol=length(Cp.v)+5,
                       dimnames=list(rownames(BIN.MX), c("mutbin", "Tmut", "Nmsite", 
                                                         "TmutDIVNmsite", "Nmsitenorm",
                                                         sort(Cp.v)))
                       )
    MUTCP.DF[[chr]] <- cbind.data.frame(chr=chr, bins=as.numeric(rownames(BIN.MX)),
                                 MUTCP.DF[[chr]], stringsAsFactors=FALSE)
    
    # Assign bin to Cp (involved in at least 1 contact for that Cp in any tissue)
    for(Cp in Cp.v){
      Cp.ind <- grep(x=colnames(BIN.MX), pattern=paste0("s_Cp_", Cp, "_ct"))
      if(length(Cp.ind)!=length(Cp.v)){ stop(paste0(chr, ": Checkpoint 1.")) }
      # Rowsums only for at least 2 rows
      MUTCP.DF[[chr]][,as.character(Cp)] <- as.numeric(rowSums(x=BIN.MX[,Cp.ind])>0)
    }
    
    chr.TF <- MUTBIN.DF$chr==chr
    bin.mut <- as.character(MUTBIN.DF$bin[chr.TF])
    MUTCP.DF[[chr]][bin.mut,c("Tmut", "Nmsite", "TmutDIVNmsite",  "Nmsitenorm")] <- MUTBIN.DF[chr.TF, c("bin", "Tmut", "Nmsite", "TmutDIVNmsite",  "Nmsitenorm")]
    MUTCP.DF[[chr]][bin.mut,"mutbin"] <- 1
   
    print(paste0(chr, " done!"), quote=FALSE)
    rm(BIN.MX); gc()
    
  } # chr.v for loop end

  MUTCP.DF <- do.call("rbind.data.frame", c(MUTCP.DF, stringsAsFactors=FALSE))
  rownames(MUTCP.DF) <- NULL
  
  mutcp.id <- paste0(gcb, "_", data.id, "_", src.id, "_", mut.id)
  save(MUTCP.DF, file=paste0(out.dir, "/", mutcp.id, "_mutCalcPerCp.RData"))
 
  #-------------------Boxplots
  x <- reshape2::melt(data=MUTCP.DF[,-(1:2)], id.vars=c("mutbin", "Tmut", "Nmsite", "TmutDIVNmsite",  "Nmsitenorm"))
  rm(MUTBIN.DF, MUTCP.DF); gc()
  # Take only bins with Cp >= 1
  x <- x[x$value==1,colnames(x)!="value"]
  colnames(x) <- c(colnames(x)[1:5], "Cp")
  
  for(calc in calc.v){
    
    bin.TF <- is.finite(x[[calc]])
    binPerCp <- table( x$Cp[bin.TF] )
    binmutPerCp <- table(x$Cp[x$mutbin==1 & bin.TF])
    percbinmut <- binmutPerCp/binPerCp*100
    binPerCp <- paste0("\n binPerCp1to21_", paste(x=binPerCp, collapse="_"))
    binmutPerCp <- paste0("\n binmutPerCp1to21_", paste(x=binmutPerCp,collapse="_"))
    percbinmut <- paste0("\n %binmutPerCp1to21_",
                         paste(x=round(percbinmut, digits=3), collapse="_"))
    
    # By default, ignore missing values in either the response or the group.
    eval(parse(text=paste0(
      'boxplot(', calc, '~Cp, outline=FALSE, data=x, xlab="Cp", ylab=calc, boxwex=0.6, 
              cex.axis=1.2, col="#FDC776", cex.main=0.3,
              main=paste0(mutcp.id, "_', calc, '_", binPerCp, binmutPerCp, 
                          percbinmut))'
    )))
    
  }
  
  print(paste0(mut, " done!"), quote=FALSE)
  
} # mut.v for loop end

dev.off()

# rm(list=ls()); gc()

