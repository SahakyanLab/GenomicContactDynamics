y <- colorRampPalette(rev(brewer.pal(11, "Spectral")))(21)
x <- col2rgb( y )
colnames(x) <- y
write.csv(x, file=paste0(output.dir, "/revRainbow_rgb.csv"), 
          quote=FALSE, row.names=TRUE )