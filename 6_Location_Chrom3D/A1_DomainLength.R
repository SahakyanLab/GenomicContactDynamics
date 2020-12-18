################################################################################
# Histogram of lengths of domains (represented as beads). Domains are TADs, LADs
# and the gaps between them. 
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
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/6_Location_Chrom3D"
    os = "Mac"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
model.id = "H1-hESC_LMNB1_hg38" #"IMR90_LMNB1_GSE49341_hg19" | "H1-hESC_LMNB1_hg38"
out.dir = paste0(wk.dir, "/out_DomainLength")
### OTHER SETTINGS #############################################################
domainfile = paste0(wk.dir, "/0_gtrack/", model.id, "/IMR90_inter_intra_chr_w_LADs.haploid.gtrack")
domainfile = paste0(wk.dir, "/0_gtrack/", model.id, "/H1-hESC_inter_intra_chr_w_LADs.haploid.gtrack")
out.name = paste0(model.id, "_Chrom3D_DomainLengths")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# load domain file
DOMAIN.coord <- fread(file=domainfile, header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)
DOMAIN.coord$domainlen <- as.numeric(DOMAIN.coord$end-DOMAIN.coord$start)
mean.val.bp <- mean(DOMAIN.coord$domainlen)
DOMAIN.coord$domainlen <- log10(DOMAIN.coord$domainlen)

if( sum(DOMAIN.coord$domainlen<0)!=0 ){
  stop("Some DomainStart > DomainEnd")
}

mean.val <- mean(DOMAIN.coord$domainlen)
ggplot(data=DOMAIN.coord, aes(x=domainlen)) +
  geom_histogram(fill="#55bde6", colour="black", binwidth=0.1) +
  geom_vline( linetype="dashed", colour="black", size=0.7, 
              aes(xintercept=log10(4e4)) ) +
  geom_vline( linetype="dashed", colour="darkred", size=0.7, 
              aes(xintercept=mean.val) 
              ) +
  scale_y_continuous(limits=c(0,700)) +
  # The warning is due to this set limit which introduces values not present
  # in dataframe, should not affect the plot
  scale_x_continuous(limits=c(4, 8)) +
  labs(title=paste0(out.name, "_", length(DOMAIN.coord$domainlen), 
                    "totaldomains_binwidth=0.1_meanLen=", mean.val.bp, "bp"),
       x=bquote(bold("log"["10"]~"("~"L"^"dom"~")")), 
       y=expression(bold( "N"["dom"] ))
       ) +
  bgr2 +
  theme(plot.title=element_text(size=10))

ggsave(filename=paste0(out.dir, "/", out.name, ".pdf"), units="in", 
       width=10, height=10)

# rm(list=ls()); gc()
