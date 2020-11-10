#pdf(file=paste0(out.dir, "/", chr, "_", gcb, "_gapdist_Cp21.pdf"), width=20, height=10)
#par(mfrow=c(2,4))

#---------------------------------------
# Chr1

# Ranges: 0-10, 20-35, >200
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylim=c(0, 50), ylab=bquote(bold( .(ylb) ))
)

# Zoom-in to 0-50, 0-10 range and 20-35
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(0,50), ylim=NULL, breaks=1:220
     
)

# Zoom-in to 0-10
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(0,10), ylim=NULL, breaks=1:220
     
)

# Small freq values at 0-10
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(0,10), ylim=c(0,10), breaks=1:220
     
)

# Zoom-in to 20-35
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(20,35), ylim=c(0,50), breaks=1:220
     
)

# Small freq values at 20-35
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(20,35), ylim=c(0,10), breaks=1:220
     
)

# Zoom-in to >200
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(150,220), ylim=c(0,5), breaks=1:220
     
)
#---------------------------------------
# Chr4

# Ranges 0-15, 40-55, 180-200
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylim=c(0, 50), ylab=bquote(bold( .(ylb) ))
)

# Zoom-in to 0-15, 40-55
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(0,50), ylim=NULL, breaks=1:220
     
)

# Zoom-in to 0-15
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(0,15), ylim=NULL, breaks=1:220
     
)

# Small freq values at 0-15
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(0,15), ylim=c(0,10), breaks=1:220
     
)

# Zoom-in to 40-55
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(40,55), ylim=c(0,10), breaks=1:220
     
)

# Zoom-in to 180-200
hist(gap, col="deepskyblue3", probability=prob , plot=TRUE,
     main=paste0(chr, "_", gcb, "_", num.ij, "ij_", "Cp=21"), 
     xlab=bquote(bold("Gap,"%*%~10^6)),
     ylab=bquote(bold( .(ylb) )),
     xlim=c(180,200), ylim=c(0,5), breaks=1:220
     
)

