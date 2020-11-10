################################################################################
# Make scatter plot of gene lengths vs. Cp
# Mac, R/3.5.2
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lencp.dir = paste0(wk.dir, "/out_geneLength")
out.dir = paste0(wk.dir, "/out_geneLength/line")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
out.id <- "length" # "length"| "%R"
config <- list(
  # c(1-colour, 2-legend.name)
  TRANSCRIPT.L = c("#FDC776", "tr"),
  INTRONS.dev.EXONS = c("#f0dab9", "int/ex"),
  EXONS.L = c("#e25d6c", "ex"),
  MEAN.EXON.L = c("#f0a5ae", "mean ex"),
  INTRONS.L = c("#3288BD", "int"),
  MEAN.INTRON.L = c("#add6f0", "mean int" ) ,
  REP.PERC.TR.L = c("#FDC776", "%R tr"),
  REP.PERC.EX.L = c( "#e25d6c", "%R ex" ),
  REP.PERC.INT.L = c( "#3288BD", "%R int" )
)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(out.id=="%R"){
  lengths.v <- c("REP.PERC.TR.L", "REP.PERC.EX.L", "REP.PERC.INT.L")
} else if(out.id=="length"){
  lengths.v <- c("TRANSCRIPT.L", "INTRONS.dev.EXONS",
                 "EXONS.L", "MEAN.EXON.L",
                 "INTRONS.L", "MEAN.INTRON.L")
} else {
  stop("Invalid out.id input.")
}

id <- paste0("hg19anno_NM_", gcb)

DF <- list()

for(len in lengths.v){
  
  # Load LENCP.DF
  # Actual length in LENCP.DF is LENCP.DF$L x 10^3
  load(file=paste0(lencp.dir, "/", id, "_", len, ".RData"))
  LENCP.DF <- aggregate(formula=L~cp, data=LENCP.DF, FUN=mean)
  # Convert to length to fold change (reference is Cp=1)
  #LENCP.DF$L <- log2(LENCP.DF$L/LENCP.DF[LENCP.DF$cp==1,"L"])
  DF[[len]] <- cbind(L.name=rep(len), LENCP.DF)
  
  print(len, quote=FALSE)
  
} # lengths.v for loop end

DF <- do.call("rbind", DF)
rownames(DF) <- NULL

config <- do.call("rbind", config)
colnames(config) <- c("colour", "legendname")

DF$L.name <- factor(DF$L.name, levels=c(lengths.v))
DF$cp <- as.numeric(DF$cp)

# Plot
ggplot(data=DF, aes(x=cp, y=L, group=L.name)) +
  geom_line(aes(colour=factor(L.name)), size=2) +
  geom_point(aes(colour=factor(L.name)), size=2.5) +
  scale_colour_manual(values=config[,"colour"], 
                      labels=config[,"legendname"]
) +
  scale_x_continuous(breaks=1:21,
                     labels=as.vector(rbind(seq(from=1, to=21, by=2),
                                            rep("")
                                            )
                                      )[-22]
                     ) +
  labs(title=id,
       x=expression(bold("c"["p"])),
       y=bquote(bold("log"["2"]~"(Mean"~.(out.id)~"FC)")),
       colour="") +
  bgr2 +
  theme(legend.text=element_text(size=15, face="bold"),
        legend.title=element_text(size=20, face="bold"))

out.id <- ifelse(out.id=="%R", "repPerc", out.id)

ggsave(filename=paste0(out.dir, "/", id, "_", out.id, "_lineplot.pdf"),
       width=10, height=10, units="in")

# rm(list=ls())