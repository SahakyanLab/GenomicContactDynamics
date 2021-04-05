################################################################################
# Plot mutation calculations vs. Cp as boxplot and scatterplot. Also calculate
# significance of values across Cps relative to Cp=1.
# Generate plot per mutation type and signature. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
# Expands warnings
options(warn=1)

whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    binmx.dir = "/Users/ltamon/DPhil/GCD_polished/7_FeaturePermutation/binmx/out_bindata_1perc_HiCNorm"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
    binmx.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/20_ChromFeatAssoc/binmx/out_bindata_1perc_HiCNorm"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
mutbin.dir = paste0(wk.dir, "/out_mutCalcPerBin_KEEP")
out.dir = paste0(wk.dir, "/out_plotsVsCp")
sigExposure.file = paste0(data.dir, "/signal_mutSig/out_samplesForSignature/donorlist_signatureExposurePERCENT.csv")
### OTHER SETTINGS #############################################################
data.id = "donor_centric_PCAWG" # "CosmicNCV", "donor_centric_PCAWG"
src.id = "Hg19" # "Hg19" | "hg38ToHg19"
bin.len = 40000L

# BIN.MX
gcb = "min2Mb"
Cp.v = 1:21

# Calculation in MUTBIN.DF
calc = "numWTSEQ" 
aggregatefunx = "mean" # String of a function name e.g. "mean", "median"

#-------------------Per plot

mut.v = "All" #c("All", "C>A", "C>G", "C>T", "T>A", "T>C", "T>G")

# Load signature exposure table to get signatures
sigExposure.df <- read.csv(file=sigExposure.file, header=T, stringsAsFactors=F)
sigExposure.df <- sigExposure.df[,-c(1,3:7)]
print(paste0("sigExposure.df has ", ncol(sigExposure.df[,-1]), " signatures..."), quote=FALSE)

SIG.v = c("RefSig.MMR1_RefSig.MMR2", 
          colnames(sigExposure.df)[colnames(sigExposure.df)!="alt.ID"])

#-------------------Combined in one plot

# Should match levels of ncv.df$location 
loc.v = c("exon", "intron", "intergenic", "intron_intergenic", 
          "exon_intron_intergenic_intergenic.excl")

sigEpLim.v = c("nosampfilter", "-1_0", "10_100", "30_100", "50_100") 
#"nosampfilter" - skip filtering based on signature exposure, use all samples in the dataset 

plotOnly = F
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(compiler)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(ggsci)
source(paste0(wk.dir, "/lib/makeMUTBINDFtoMUTCPDFperMUT.R"))
source(paste0(wk.dir, "/lib/aggregateDF.R"))
source(paste0(wk.dir, "/lib/identifyAltHyp.R"))
source(paste0(wk.dir, "/lib/doMannWhitney.R"))
source(paste0(wk.dir, "/lib/makebp.R"))
source(paste0(wk.dir, "/lib/makeScatter.R"))
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Generate combinations of mut, sig.v, sigEpLim.v, loc.v
comb.mx <- data.matrix(expand.grid(mut=1, SIG=1:length(SIG.v), sigEpLim=1:length(sigEpLim.v), 
                                   loc=1:length(loc.v), stringsAsFactors=F))

# Generate combinations of mut, sig.v
combms.mx <- data.matrix(expand.grid(mut=1, SIG=1:length(SIG.v), stringsAsFactors=F))
ms.len <- nrow(combms.mx)

for(mut in mut.v){
  
  mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=T)
  out.id <- paste0(gcb, "_", data.id, "_", src.id, "_", mut.id, "_", calc, "_", aggregatefunx)
  
  p.lst   <- list()
  pfc.lst <- list()
  
  if(plotOnly==F){
    
    pdf(file=paste0(out.dir, "/", out.id, "_boxplots.pdf"), height=25, width=20)
    par(mfrow=c(5,4))
    
    P.DF <- list()
    PFC.DF <- list()
    
  } else {
    
    # P.DF
    load(file=paste0(out.dir, "/", out.id, "_scatplot.RData"))
    # PFC.DF
    load(file=paste0(out.dir, "/", out.id, "_FC_scatplot.RData"))
    
  }
    
  #ms.len <- 2 # REMOVE
  for(ms in 1:ms.len){
    
    ms.v <- as.numeric(combms.mx[ms,])
    combsl.mx <- comb.mx[ comb.mx[,"mut"]==ms.v[1] & comb.mx[,"SIG"]==ms.v[2], 
                          c("sigEpLim", "loc"), drop=F ]
    #mut <- mut[ ms.v[1] ]
    #mut.id <- gsub(x=mut, pattern=">", replacement="To", fixed=T)
    SIG <- SIG.v[ ms.v[2] ]
    
    ms.id <- paste0(gcb, "_", data.id, "_", src.id, "_", mut.id, "_", SIG)
    
    print(paste0(ms.id, "..."), quote=F)
    
    if(plotOnly==F){
      
      DF     <- list()
      DF.fc  <- list()
      sl.len <- nrow(combsl.mx)
      #sl.len <- 2 # REMOVE
      for(sl in 1:sl.len){
        
        sl.v <- as.numeric(combsl.mx[sl,])
        sigEpLim <- sigEpLim.v[ sl.v[1] ]
        loc <- loc.v[ sl.v[2] ]
        
        mutcp.id <- paste0(data.id, "_", src.id, "_", mut.id, "_", SIG, "_", loc, 
                           "_sigEperclimits_", sigEpLim)
        mutcp.file <- paste0(mutbin.dir, "/", mutcp.id, "_mutCalcPerBin.RData")
        if( file.exists(mutcp.file) ){
          load(file=mutcp.file)
        } else {
          next
        }
        
        # Generate MUTCP.DF
        MUTCP.DF <- makeMUTBINDFtoMUTCPDFperMUT(binmx.dir=binmx.dir, 
                                                gcb=gcb, Cp.v=Cp.v, 
                                                bin.len=bin.len,
                                                MUTBIN.DF=MUTBIN.DF)
        rm(MUTBIN.DF); gc()
        
        # Collect data for mut-sig plot
        MUTCP.DF <- MUTCP.DF[,c("mutbin", calc, as.character(Cp.v))]
        x <- reshape2::melt(data=MUTCP.DF, id.vars=c("mutbin", calc), stringsAsFactors=F)
        rm(MUTCP.DF); gc()
        
        colnames(x)[colnames(x)=="variable"] <- "ind"
        
        # Take only bins with Cp >= 1
        x <- x[x$value==1,colnames(x)!="value"]
        
        if( !identical(levels(x$ind), as.character(Cp.v) ) ){
          stop(paste0(mutcp.id, ": Cp levels wrong." ))
        }
        
        # Make boxplot
        BINMUT <- makebp(df=x, calc=calc, xlab="Cp", ylab=calc, addjitter=F, plot.id=mutcp.id)
        
        # x$ind left as factor for aggregateDF()
        DF[[mutcp.id]] <- aggregateDF(ind=x$ind, values=x[[calc]], FUN=aggregatefunx)
        DF[[mutcp.id]] <- cbind.data.frame(mut.id=rep(mut.id), SIG.id=rep(SIG), loc.id=rep(loc), 
                                           sigEpLim.id=rep(sigEpLim), DF[[mutcp.id]], 
                                           percbinmut=as.numeric(BINMUT$percbinmut[ DF[[mutcp.id]]$ind ]),
                                           binPerCp=as.numeric(BINMUT$binPerCp[ DF[[mutcp.id]]$ind ]),
                                           stringsAsFactors=F)
        rm(BINMUT)
        
        # Mann-Whitney
        x$ind <- as.character(x$ind)
        DF[[mutcp.id]] <- doMannWhitney(x=x, df=DF[[mutcp.id]], calc=calc)
        rownames(DF[[mutcp.id]]) <- DF[[mutcp.id]]$ind
        
        # Fold change version relative to Cp=1
        DF.fc[[mutcp.id]] <- DF[[mutcp.id]]
        DF.fc[[mutcp.id]]$values <- log2(DF[[mutcp.id]]$values/DF[[mutcp.id]]["1","values"])
        
        print(paste0(sigEpLim, " ", loc, " done!"), quote=F)
        rm(sigEpLim, loc, mutcp.id, sl.v, x); gc()
        
      } # sl.len for loop end
      
      DF <- do.call("rbind.data.frame", c(DF, stringsAsFactors=F))
      DF.fc <- do.call("rbind.data.frame", c(DF.fc, stringsAsFactors=F))
      rownames(DF) <- rownames(DF.fc) <- NULL
      
      P.DF[[ms.id]] <- DF
      PFC.DF[[ms.id]] <- DF.fc
      
      rm(sl.len)
      
    } else {
      
      DF <- P.DF[ P.DF==mut.id & P.DF$SIG.id==SIG, ]
      DF.fc <- PFC.DF[ PFC.DF==mut.id & PFC.DF$SIG.id==SIG, ]
      
    }
    
    #-------------------PLOT
    
    p.lst[[ms.id]] <- makeScatter(df=DF, calc=calc, yint=NULL, plot.id=ms.id)
    pfc.lst[[ms.id]] <- makeScatter(df=DF.fc, calc=calc, yint=0, plot.id=paste0(ms.id, "_FCreltoCp=1"))
    
    rm(SIG, ms.v, combsl.mx, ms.id, DF, DF.fc); gc()
    
  } # ms.len for loop end
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=5, ncol=4, legend=NULL)
  ggexport(p.arr, height=50, width=60, filename=paste0(out.dir, "/", out.id, "_scatplot.pdf"))
  
  p.arr <- ggarrange(plotlist=pfc.lst, nrow=5, ncol=4, legend=NULL)
  ggexport(p.arr, height=50, width=60, filename=paste0(out.dir, "/", out.id, "_FC_scatplot.pdf"))
  
  if(plotOnly==F){
    
    dev.off()
    
    P.DF <- do.call("rbind.data.frame", c(P.DF, stringsAsFactors=F))
    PFC.DF <- do.call("rbind.data.frame", c(PFC.DF, stringsAsFactors=F))
    rownames(P.DF) <- rownames(PFC.DF) <- NULL
    
    if( any(!is.character(c(P.DF$ind, PFC.DF$ind))) ){
      stop(paste0(out.id, ": Attribute of indices not character."))
    }
    P.DF$ind <- as.numeric(P.DF$ind)
    PFC.DF$ind <- as.numeric(PFC.DF$ind)
    
    save(P.DF, file=paste0(out.dir, "/", out.id, "_scatplot.RData"))
    save(PFC.DF, file=paste0(out.dir, "/", out.id, "_FC_scatplot.RData"))
    
  }

  print(paste0(mut.id, " done!"), quote=F)
  
  rm(p.arr, p.lst, pfc.lst, P.DF, PFC.DF, mut, mut.id, out.id)
  
} # mut.v for loop end

# rm(list=ls()); gc()