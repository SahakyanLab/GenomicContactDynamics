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
  HYB.MX = HYB.MX,
  # Corresponding order to HYB.MX (column)
  cp.v = PERSIST.MX$ntis,
  # Corresponding order to HYB.MX (row)
  label.v = c("T.H.E. Gfree, kcal/mol", 
             "T.H.E. sd(discordance)", 
             "T.H.E. sum(abs(discordance))")
){
  
  feat.v <- dimnames(HYB.MX)[[1]]
  feat.v.len <- length(feat.v)
  names(label.v) <- feat.v
  
  # Numbers
  tot.ij <- length(HYB.MX[1,])
  incl.ij.TF <- !is.na(HYB.MX[1,])
  incl.ijPerc <- round(x=sum(incl.ij.TF)/tot.ij*100, digits=2)
  percPerCp <- table(cp.v[incl.ij.TF])/tot.ij*100
  percPerCp <- round(x=percPerCp[as.character(1:21)], digits=4)
  
  for(i in 1:feat.v.len){
    
    par(oma=c(0,0,0,0), mar=c(5, 5, 4, 2) + 0.1, mgp=c(3, 1.2, 0))
    pdf(file=paste0(out.dir, "/", out.name, "_", feat.v[i], "_indivBP.pdf"),
        width=10, height=8)
    
    df <- cbind(value=HYB.MX[i,], cp=cp.v)
    
    boxplot(value~cp, outline=FALSE, data=df, xlab="", ylab="", boxwex=0.6,
            cex.axis=1.2, col="#FDC776", xaxt="n")
    rm(df); gc()
    # X-axis
    axis(side=1, at=1:21, labels=1:21, cex.axis=1.2)
    # X-axis title
    mtext(side=1, text=expression(bold( "c"["p"] )), line=3.2, cex=2)
    # Y-axis title
    mtext(side=2, text=bquote(bold( .(label.v[i]) )), line=2, cex=2)
    # Plot title
    mtext(side=3, text=paste0(out.name, "_", tot.ij, "totij_", incl.ijPerc, "%nonNA"), line=2, cex=1)
    mtext(side=3, text=paste(percPerCp, collapse=" "), line=0.2, cex=0.5)
    
    dev.off()
    
  }
  
}

