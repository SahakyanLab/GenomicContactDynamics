################################################################################
# Generate plots showing kmer-based discordance score of HiC contacts for each
# chromosome
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
# library(RColorBrewer)
### FUNCTION ###################################################################
HiCHybridPlot <- function(
  out.dir = "/dir",
  out.name = paste0(chr, "_" , gcb, "_kmer", k, affix),
  HYB.MX = t(CII.MX),
  # Corresponding order to HYB.MX (column)
  cp.v = PERSIST.MX$ntis,
  # Corresponding order to HYB.MX (row)
  label.list = list(bquote(bold( "Gfree, kcal/mol" )),
                    bquote(bold( "s ("~"C"["||"]~")" )),
                    bquote(bold( "C"["||"] ))
                    )
){
  
  
  tot.ij <- length(HYB.MX[1,])
  nonNA.ij.TF <- !is.na(HYB.MX[1,])
  percPerCp <- table(cp.v[nonNA.ij.TF])/tot.ij*100
  percPerCp <- format(x=percPerCp[as.character(1:21)], digits=2, scientific=TRUE)
  percPerCp <- paste(percPerCp, collapse=" ")
  
  plottitle <- paste0(out.name, "_", format(x=as.numeric(tot.ij), digits=3,
                                            scientific=TRUE), 
                      "totij_", format(x=100*sum(nonNA.ij.TF)/tot.ij, 
                                       scientific=TRUE, digits=2), 
                      "%nonNA\n", percPerCp)
  

  feat.v <- dimnames(HYB.MX)[[1]]
  feat.v.len <- length(feat.v)
  label.list <- as.vector(label.list)
  names(label.list) <- feat.v
  
  for(i in 1:feat.v.len){
    
    par(oma=c(0,0,0,0), mar=c(5, 5, 4, 2) + 0.5, mgp=c(3, 1.2, 0))
    prefix <- strsplit(x=out.name, split="\\n")[[1]][1]
    pdf(file=paste0(out.dir, "/", prefix, "_", feat.v[i], "_indivBP.pdf"),
        width=10, height=8)
    
    df <- cbind(value=HYB.MX[i,], cp=cp.v)
    
    boxplot(value~cp, outline=FALSE, data=df, xlab="", ylab="", boxwex=0.6,
            cex.axis=1.2, col="#FDC776", xaxt="n")
    rm(df); gc()
    # X-axis
    axis(side=1, at=1:21, labels=1:21, cex.axis=1.2)
    # X-axis title
    mtext(side=1, text=expression(bold( "c"["p"] )), line=3.1, cex=2)
    # Y-axis title
    mtext(side=2, text=label.list[[i]], line=2.0, cex=2)
    # Plot title
    mtext(side=3, text=plottitle, line=2, cex=0.7)
    
    dev.off()
    
  }
  
}

