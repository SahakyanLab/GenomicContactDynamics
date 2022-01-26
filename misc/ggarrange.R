myplots.arranged <- ggpubr::ggarrange(plotlist=myplots.list, nrow=1, ncol=3)
ggpubr::ggexport(myplots.arranged, width=30, height=10,
                 filename=paste0(output.dir, "/bp_", genome.ver, "_min", gcb, 
                                 "Mb_chrALL_TxSizeVsPersist.pdf"))