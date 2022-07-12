################################################################################
# Generate bootstrapped estimates for mean, sd, var, sd/mean, fraction of NE, 
# LE, ME, HE of gene expression expected by random chance. Bootstrapping is done
# per hub and per tissue because n or the number of genes with data per hub per
# tissue varies. A random sample of n size is drawn from the genes of the 
# corresponding hub chr only considering genes with expression data for the
# corresponding tissue. Also, To simplify the process of drawing the random 
# sample, I only considered genes whose transcripts are located in the one chr.
# Values in the random sample below the set expression cut-off
# (in our case 0.5) is converted to 0. Bootstrapping is always done 10000 times. 
# Using the seed0, 10000 seeds are generated, each will be used for each iteration.
# In other words, the same set of 10000 integers are used in the bootstrappping.
# However, values of the random samples will be different because n varies across
# hubs and tissues as well as the set of genes that have expression data. 
# When SDdivMEAN NaN (when MEAN=SD=0), assign to 0 because they should not be 
# excluded since they are data. 0 because when all genes are not expressed then 
# there is no variation.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
# Data on genes forming hubs and the n per hub
exprData.dir = paste0(wk.dir, "/out_cleanExprData")
exprData.suffix = "cutoff0_LTr_ALL"
hub.dir = paste0(wk.dir, "/out_hubfile/All_LTr")
hubinfo.dir = paste0(wk.dir, "/out_hubinfo/All_LTr")
out.dir = paste0(wk.dir, "/out_nonHubBoot/All_LTr/exprCutoff_05")
### OTHER SETTINGS #############################################################
src.id = "data2"
gcb = "min2Mb"
hub.id = "All_topCP3_gapBin50"
expr.cutoff = 0.5

# BOOTSTRAP parameters
iter = 10000
nCPU = 3
seed0 = 763
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
exprDataPath <- paste0(exprData.dir, "/expr_", src.id, "_", exprData.suffix, ".csv")
nhub.df <- read.csv(file=paste0(hubinfo.dir, "/", gcb, "_", hub.id, "_", src.id, 
                                "_nvalues.csv"), stringsAsFactors=F) 

HUB <- list.files(path=hub.dir, pattern=hub.id, recursive=F)
HUB <- HUB[grepl(x=HUB, pattern=gcb)]
HUB <- gsub(x=HUB, pattern=".csv", replacement="", fixed=T)
HUB.len <- length(HUB)
print(paste0(HUB.len, " hubs."), quote=F)

expr.df <- read.csv(file=exprDataPath, header=T, stringsAsFactors=F)
expr.df <- expr.df[!is.na(expr.df$chr),]
tiss.v <- colnames(expr.df[,!colnames(expr.df)%in%c("chr", "Gene.Name")])
tiss.v.len <- length(tiss.v)

stat.v <- c("MEAN", "sd.MEAN", "MEDIAN", "sd.MEDIAN", "SD", "sd.SD", 
            "VAR", "sd.VAR", "NE", "sd.NE", "LE", "sd.LE", "ME", 
            "sd.ME", "HE", "sd.HE", "SDdivMEAN", "sd.SDdivMEAN")
stat.v.len <- length(stat.v)

set.seed(seed0)
# Generate seeds per hub
seed.hub <- sample(x=1:1000, size=HUB.len, replace=F)

for(h in 1:HUB.len){
  
  set.seed(seed.hub[h])
  # Generate seeds per hub, per tissue
  seed.tiss <- sample(x=1:1000, size=tiss.v.len, replace=F)
  names(seed.tiss) <- tiss.v
  
  hub <- HUB[h]
  hub.ind <- which(nhub.df$hub==hub) 
  chr <- nhub.df[hub.ind, "chr"]
  BOOT.MX <- matrix(data=NA, nrow=tiss.v.len, ncol=stat.v.len+1,
                    dimnames=list(tiss.v, c("n", stat.v)))
  BOOT.MX[,"n"] <- 0
  
  for(tiss in tiss.v){
    
    set.seed( unname(seed.tiss[tiss]) )
    # Generate seeds per hub, per tissue, per iteration
    SEED.v <- sample(x=1:20000, size=iter, replace=F)
    
    n <- nhub.df[hub.ind, tiss]
    # Just diregard genes with transcripts on different chr, there are 
    # enough genes per chr anyway
    val <- expr.df[expr.df$chr==paste0(chr, "."),tiss]
    val <- val[!is.na(val)]
    
    if(n==0){
      
      print(paste0(n, " ", tiss, " skipped!"), quote=F)
      next
      
    } else {
      
      toExport <- c("SEED.v", "val", "n")
      #### PARALLEL EXECUTION #########
      boot <- foreach(itr=isplitVector(1:iter, chunks=nCPU), 
                      .inorder=T, .combine="rbind",
                      .export=toExport, .noexport=ls()[!ls()%in%toExport]
                      
      ) %op% {
        
        chunk <- sapply(X=itr, simplify=F, FUN=function(i){
          
          set.seed(SEED.v[i])
          samp <- sample(x=val, size=n, replace=T)
          samp[samp<expr.cutoff] <- 0
          samp <- c(MEAN=mean(samp), MEDIAN=median(samp), SD=sd(samp), VAR=var(samp),
                    # Based on https://www.ebi.ac.uk/gxa/help/index.html
                    NE=sum(samp==0)/n, 
                    LE=sum(samp>=expr.cutoff & samp<=10)/n,
                    ME=sum(samp>10 & samp<=1000)/n,
                    HE=sum(samp>1000)/n)
          #if( sum(samp[4:7])!=1 ){ stop("Sum not 1.") }
          samp <- c(samp, SDdivMEAN=as.numeric(samp["SD"]/samp["MEAN"]))
          return(samp)
          
        })
        
        return(do.call("rbind", chunk))
        
      }
      rm(val)
      ### END OF PARALLEL EXECUTION ###
      
      # SDdivMEAN NaN (when MEAN=SD=0), assign to 0 because they
      # should not be excluded since they are data. 0 because
      # when all genes are not expressed then there is no variation.
      boot[ is.nan(boot[,"SDdivMEAN"]),"SDdivMEAN" ] <- 0
      
      temp <- apply(X=boot, MARGIN=2, FUN=function(stat){
        c(mean(stat), sd(stat))
      })
      
      rm(boot)
      
    }
    
    BOOT.MX[tiss,] <- c(n, as.vector(temp))
    rm(temp); gc()
    print(paste0(n, " ", tiss, " done!"), quote=F)
    
    rm(SEED.v)
    
  } # tiss.v for loop end
  
  rm(seed.tiss)
  
  if( any(is.nan(BOOT.MX)) ){
    stop(paste0(hub, ": NaN values in BOOT.MX"))
  }
  save(BOOT.MX, file=paste0(out.dir, "/", hub, "_", src.id, 
                            "_seed", seed0, "_iter", iter, ".RData"))
  print(paste0(hub, " done!"), quote=F)
  
} # HUB.len for loop end

# rm(list=ls()); gc()



