lib = "~/DPhil/lib"

library(ggplot2)
library(RColorBrewer)
source(paste0(lib, "/GG_bgr.R"))

df <- data.frame(Cp=1:21, x=1:21, y=1)


cols <- colorRampPalette(rev(brewer.pal(n=11, name="Spectral")))(21)

p <- ggplot(data=df, aes(x=x, y=y)) +
  geom_point(aes(col=Cp)) + 
  scale_colour_gradientn(colours=cols) +
  bgr2

ggsave(file=paste0("./rainbow_legend_continuous.pdf"), plot=p, width=10, height=10)

