################################################################################
# Function to generate CORCP.DF
### FUNCTION ###################################################################
generateCORCPDF <- function(chr.v, 
                            pairCor.dir, src.id, expr.cutoff, gene.id, cor.meth, 
                            pairCp.dir, gcb, bin.len, 
                            pairHub.dir, hub.id,
                            out.dir, out.name
                            ){
  
  CORCP.DF <- list()
  for(chr in chr.v){
    
    # PAIRCOR.MX
    load(paste0(pairCor.dir, "/expr_", src.id, "_cutoff", expr.cutoff, "_", gene.id, 
                "_", chr, "_", cor.meth, "_coexpression.RData"))
    
    PAIRCOR.MX <- PAIRCOR.MX[order(PAIRCOR.MX[,1], PAIRCOR.MX[,2]),]
    
    load(paste0(pairCp.dir, "/", chr, "_", gcb, "_", gene.id, "_binlen", bin.len, 
                "_pairMaxCp_coexpression.RData"))
    dimnames(PAIRCP.MX)[[2]] <- c("g1.ind", "g2.ind", "maxCp")
    
    # PAIRCP.MX
    PAIRCP.MX <- PAIRCP.MX[order(PAIRCP.MX[,1], PAIRCP.MX[,2]),]
    
    fle <- paste0(pairHub.dir, "/", chr, "_HUBID", hub.id, "_", gene.id, 
                  "_pairHub_coexpression.RData")
    if( file.exists(fle) ){
      load(fle)
    } else {
      PAIRHUB.MX <- cbind(PAIRCOR.MX[,1:2], WithinHub.TF=NA)
    }
    
    dimnames(PAIRHUB.MX)[[2]] <- c("g1.ind", "g2.ind", "hubTF")
    
    # PAIRHUB.MX
    PAIRHUB.MX <- PAIRHUB.MX[order(PAIRHUB.MX[,1], PAIRHUB.MX[,2]),]
    
    check1.TF <- !identical( PAIRCOR.MX[,1:2], PAIRCP.MX[,1:2] )
    check2.TF <- !identical( PAIRCOR.MX[,1:2], PAIRHUB.MX[,1:2] )
    
    #
    if(check1.TF | check2.TF){
      stop("Gene pair indices of PAIRCOR.mx and PAIRCP.MX not identical.")
    } else {
      
      CORCP.DF[[chr]] <- cbind(PAIRCOR.MX[,c("g1.ind", "g2.ind")],
                               corval=PAIRCOR.MX[,"cor"], 
                               maxCp=PAIRCP.MX[,"maxCp"],
                               hubTF=PAIRHUB.MX[,"hubTF"],
                               frTisspairAtleast1NA=PAIRCOR.MX[,"frTisspairAtleast1NA"])
    }
    
    rm(PAIRCOR.MX, PAIRCP.MX)
    print(paste0(chr, " data obtained!"), quote=F)
    
  }
  
  gc()
  
  CORCP.DF <- do.call("rbind.data.frame", CORCP.DF)
  rownames(CORCP.DF) <- NULL
  
  write.csv(x=stack(table(CORCP.DF$maxCp, useNA="always")), row.names=F, 
            file=paste0(out.dir, "/", out.name, "_NdatapointsPerCp.csv"))
  
  save(CORCP.DF, file=paste0(out.dir, "/", out.name, "_plot.RData"))
  
  return(CORCP.DF)
  
}

# rm(list=ls()); gc()