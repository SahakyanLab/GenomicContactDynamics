################################################################################
# Apply markov clustering on Cs and Cp-based contact matrices.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_MCL"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_MCL"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
#persist.dir = paste0(wk.dir, "/out_basePersist")
out.dir = paste0(wk.dir, "/out_mcl")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
topCP = 21; cp.v = 1:21
ct = "FC"; ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
                    "AG", "Ov", "Bl", "MesC","MSC", "NPC", "TLC", "ESC", "FC", "LC")
gap.bin = 0
e.val = 2 # MCL - Expansion
gran.val.v = arr2.repl #c(2,3,4,5,10,20,30,40,50,60,80) # MCL - Inflation
weight = "Cp" # Value | "Cp" | "Cs"
plotOnly = FALSE
makeADJMX = FALSE
combinePlot = TRUE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(expm)
library(sna)
source(paste0(lib, "/MyMCL.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
ct.id <- NULL
if( !is.null(ct) ){ ct.id <- ct; ct.id <- paste0(ct.id[!is.na(ct.id)], "_") }
out.id <- paste0(ct.id, "topCP", topCP, "_gapBin", gap.bin, "_weight", weight)

for(gran.val in gran.val.v){
  
  mcl.id <- paste0("gran", gran.val, "_expan", e.val)
  print(out.id, quote=FALSE)
  print(mcl.id, quote=FALSE)
  
  mclout.dir <- paste0(out.dir, "/", out.id, "_", mcl.id)
  if( !dir.exists(path=mclout.dir) ){
    dir.create(path=mclout.dir)
  }
  
  if(combinePlot){
    pdf(file=paste0(out.dir, "/", gcb, "_", out.id, "_", mcl.id, ".pdf"), 
        height=40, width=60)
    par(mfrow=c(4,6))
  }
  
  for(chr in chr.v){
    
    out.name <- paste0(gcb, "_", chr, "_", out.id, "_", mcl.id)
    
    if(plotOnly==FALSE){
      
      if(makeADJMX){
        load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
        #load(file=paste0(persist.dir, "/", chr, "_Persist_min2Mb_topCP4.RData"))
        if(topCP<0){ cp.v <- rev(cp.v) }
        cp.TF <- PERSIST.MX$ntis%in%tail( x=cp.v, n=abs(topCP) ) 
        gap.TF <- ((PERSIST.MX$hits$j-PERSIST.MX$hits$i)-1) >= gap.bin
        if( is.null(ct) ){
          ij.TF <- gap.TF & cp.TF
        } else if(ct%in%ct.v){
          ij.TF <- gap.TF & cp.TF & PERSIST.MX$hits[[ct]]>0
        }; rm(gap.TF, cp.TF)
        ij.df <- cbind(PERSIST.MX$hits[ij.TF,c("i", "j")], Cp=PERSIST.MX$ntis[ij.TF])
        ij.len <- nrow(ij.df)
        if(ij.len==0){
          print(paste0(chr, ": No contacts for given gap, Cp and ct (if given)."), 
                quote=FALSE)
          next
        }
        
        if(weight=="Cp"){
          ## Scale Cp values using the maximum
          #weight.v <- PERSIST.MX$ntis[ij.TF]/max(cp.v)
          weight.v <- PERSIST.MX$ntis[ij.TF]
        } else if( weight=="Cs" & !is.null(ct) ){
          weight.v <- PERSIST.MX$hits[ij.TF, ct]
          #weight.v <- weight.v/max(weight.v)
        } else if( is.numeric(weight) & weight>0 ){
          weight.v <- rep(x=weight, times=ij.len)
        } else {
          stop("Specify weight and ct accordingly.")
        }
        rm(PERSIST.MX, ij.TF); gc()
        
        ubins <- unique( c(unique(ij.df$i), unique(ij.df$j)) )
        ubins <- sort(ubins, decreasing=FALSE)
        ubins.len <- length(ubins)
        ADJMX <- matrix(data=0, nrow=ubins.len, ncol=ubins.len, 
                        dimnames=list(ubins,ubins))
        for(n in 1:ij.len){
          ij <- as.character(ij.df[n,])
          w.val <- weight.v[n]
          ADJMX[ij[1],ij[2]] <- w.val
          ADJMX[ij[2],ij[1]] <- w.val
        }
        rm(weight.v, ij, w.val, ubins); gc()
        print("ADJMX done!", quote=FALSE)
        save(ADJMX, file=paste0(out.dir, "/", gcb, "_", chr, "_", out.id, "_adjmx.RData"))
      } else {
        # LOAD ADJMX
        load(file=paste0(out.dir, "/", gcb, "_", chr, "_", out.id, "_adjmx.RData"))
      }
      print("Apply MCL...", quote=FALSE)
      MCLOUT <- MyMCL(adjMX=ADJMX, e=e.val, gran=gran.val, addLoops=TRUE, 
                      max.iter=1000, allowVC=FALSE)
      save(MCLOUT, file=paste0(mclout.dir, "/", out.name, "_mcl.RData"))
      
    } else {
      flenme <- paste0(mclout.dir, "/", out.name, "_mcl.RData")
      if(file.exists(flenme)){
        load(file=flenme)
      } else {
        print(paste0(chr, ": No MCL output for given gap and granularity."), quote=FALSE)
        next
      }
    }
    
    if(!combinePlot){
      pdf(file=paste0(mclout.dir, "/", out.name, ".pdf"), height=10, width=10)
    }
    set.seed(100) # To replicate the position of nodes 
    POS.MX <- gplot(MCLOUT$equilibrium.state, gmode="graph", displaylabels=TRUE,
                    label.cex=0.5, vertex.cex=0.5, main=out.name, cex=0.5)
    save(POS.MX, file=paste0(mclout.dir, "/", out.name, "_posmx.RData"))
    if(!combinePlot){
      dev.off()
    }
    print(paste0(chr, " done!"), quote=FALSE)
    
  }
  if(combinePlot){
    dev.off()
  }
  
} # gran.val.v for loop end

# rm(list=ls())



