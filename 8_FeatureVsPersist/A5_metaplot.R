################################################################################
# Two plots:
# a. Get mean %bin for each 40kb-bin then calculate fold change relative to mean of 
# actual HiC bin. Plot is log2(FC) X midpoint of each bin(x10^4)
# b. See what feature go down or up with increasing Cp. Per feature, plot
# mean %bin of actual HiC bin Vs contact persistence. 
# deva, R/3.6.0-newgcc, gcc/4.9.2

# ERROR: missing points in plot, this is due to preset y-axis min and max limits 
# in plotParam (i.e. -1, 1), change if necessary
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/8_FeatureVsPersist"
    data.dir = "/Users/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# Chromatin features directory
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/raw")
foifile = paste0(wk.dir, "/foifile/FINAL_foifile_priority_nperm10000_seed662_mxmskfr0_Cp21_pvalcutoff0.05_numOlapA_comOlap_edited_consistentAcrossCT")
foifile = paste0(wk.dir, "/foifile/foifile_genes")
# File of feature grouping
featgrpfile = paste0(wk.dir, "/features_group")
fetacp.dir = paste0(wk.dir, "/out_FETACP")
out.dir = paste0(wk.dir, "/out_metaplot")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
description = NULL
out.id = "genes" #"priority_Cp21_numOlapAANDcomOlap_acrossTissues"
plotOnly = FALSE
# What bin to plot for FC vs. bin position
cp = 21
# What value to use for ordering?
FCBIN.order = 1 # pos 
FCCP.order = 21 # Cp

colourBy = "foi" # "foi" | "group"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(reshape)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)
library(yarrr); base.pal <- yarrr::piratepal("basel")
source(paste0(lib, "/GG_bgr.R"))
source(paste0(wk.dir, "/lib/finaliseFOI.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# File with feature groupings
featgrp <- readLines(con=featgrpfile)
header.ind <- grep(x=featgrp, pattern=">")
# Remove ">" from the headers
featgrp[header.ind] <- gsub(x=featgrp[header.ind], pattern=">",
                            replacement="")
header.ind <- c(header.ind, length(featgrp)+1)
header.ind.len <- length(header.ind)
featgrp.v <- sapply(X=1:(header.ind.len-1), simplify=TRUE, FUN=function(x){
  v <- featgrp[(header.ind[x]+1):(header.ind[x+1]-1)]
  v <- paste("_", v, "_", sep="")
  #return( paste(v, collapse="|") )
})
names(featgrp.v) <- featgrp[header.ind[-header.ind.len]]
feat.df <- stack(featgrp.v)
feat.df$ind <- as.character(feat.df$ind)
#-------------------------------------------------------------------------------
id <- paste0("chrALL_", gcb, "_desc_", out.id)

if(plotOnly==FALSE){
  
  # List of features
  foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)
  
  # Filter feature list by description
  foi.v <- foi.v[ grepl(x=foi.v, pattern=ifelse(is.null(description), "", description), 
                        ignore.case=FALSE) ]
  foi.v.len <- length(foi.v)
  
  # Order plots by fold change certain FC value (mentioned below); ± sign is 
  # accounted in the order so that both extremes can be seen
  ORDER <- list()
  # Order plots max fold change (relative to contacting bin/ bin 0)
  ORDER$FCVSBIN <- rep(NA, times=foi.v.len)
  # Order plots by fold change (relative to Cp=1) at Cp=21
  ORDER$FCVSCP <- ORDER$FCVSBIN
  
  # FC for Cp=21 across region spanning contacting bin
  FCVSBIN <- list()
  # FC at contacting bin across Cp
  FCVSCP <- list()
  
  # Store group of each foi
  for(f in 1:foi.v.len){
    
    foi.v[f] <- tail(x=strsplit(x=foi.v[f], split="\\/")[[1]], n=1)
    foipat <- strsplit(x=foi.v[f], split="foi\\_|\\_desc")[[1]][2]
   
    # Load FETACP.MX
    load(file=paste0(fetacp.dir, "/chrALL_", gcb, "_", 
                     gsub(x=foi.v[f], 
                          pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement=""),
                     "_fetacp.RData"))
    
    # Foi name for output
    #foi.v[f] <- strsplit(x=foi.v[f], split="foi\\_|\\_desc\\_|\\.bed")[[1]][2]
    fgrp <- strsplit(x=foi.v[f], split="foi\\_|\\_desc\\_|\\.bed")[[1]][2]
    
    grp.TF <- unlist(sapply(X=feat.df$values, simplify=TRUE, FUN=function(pat){
      grepl(pattern=pat, x=paste0("_", foipat, "_"))
    }))
    if(foipat=="A_Phased_Repeat"){
      grp.TF[!is.na(grp.TF)] <- FALSE 
      grp.TF["_A_Phased_Repeat_"] <- TRUE
    }
    if(sum(grp.TF)>1){ stop( paste0(fgrp, ":Promiscous group for feature.")) }
    grp <- ifelse(sum(grp.TF)==0, "Rest", feat.df$ind[grp.TF])
    
    if(f==1){
      pos.v <- as.numeric(colnames(FETACP.MX))
      labs <- pos.v
    }
    
    # Bin FETACP.MX data (40 kb)
    df.orig <- melt.array(FETACP.MX)
    colnames(df.orig) <- c("Cp", "bin", "percbin")
   
    #---------------------------------------------------
    # FC for Cp=21 across region spanning contacting bin
    #---------------------------------------------------
    df <- df.orig[df.orig$Cp==cp,]
    df <- aggregate(formula=percbin~bin, data=df, FUN=mean)
    # Log fold change reference is value at contacting bin, pos=0
    meanCpref <- df[df$bin==0, "percbin"]
    # Fold change
    df$percbinDIVref <- df$percbin/meanCpref
    df <- cbind(foi=rep(foi.v[f]), df)
    # Value for ordering plots
    ORDER$FCVSBIN[f] <- log2(df[df$bin==FCBIN.order, "percbinDIVref"])
    FCVSBIN[[ foi.v[f] ]] <- df
    
    #---------------------------------------------------
    # FC at contacting bin across Cp
    #---------------------------------------------------
    df <- df.orig[df.orig$bin==0,]
    df <- aggregate(formula=percbin~Cp, data=df, FUN=mean)
    # Log fold change reference is Cp=1 value
    meanCp1 <- df[df$Cp==1, "percbin"]
    df$percbinDIVref <- df$percbin/meanCp1
    df <- cbind(foi=rep(foi.v[f]), df)
    # Value for ordering plots
    ORDER$FCVSCP[f] <- log2(df[df$Cp==FCCP.order, "percbinDIVref"])
    FCVSCP[[ foi.v[f] ]] <- df
  
    # Add group
    FCVSBIN[[ foi.v[f] ]] <- cbind(group=rep(grp), FCVSBIN[[ foi.v[f] ]])
    FCVSCP[[ foi.v[f] ]] <- cbind(group=rep(grp),  FCVSCP[[ foi.v[f] ]])
    
    print(paste0(foi.v[f], " data obtained!"), quote=FALSE)
    rm(FETACP.MX, df, df.orig, meanCpref, meanCp1); gc()
    
  } # foi.v for loop end
  
  METAPLOT <- list(FCVSBIN=do.call("rbind", FCVSBIN), 
                   FCVSCP=do.call("rbind", FCVSCP)) 
  # Order foi
  for(x in names(METAPLOT)){
    rownames(METAPLOT[[x]]) <- NULL
    foi.ordrd <- order(ORDER[[x]], decreasing=TRUE)
    METAPLOT[[x]][["foi"]] <- factor(x=METAPLOT[[x]][["foi"]],
                                     levels=unique(METAPLOT[[x]][["foi"]])[foi.ordrd])
  }
  save(METAPLOT, file=paste0(out.dir, "/", id, "_metaPcomb.RData"))
  
  rm(foi.ordrd); gc()
  
} else {
  # Load METAPLOT list
  load(file=paste0(out.dir, "/", id, "_metaPcomb.RData"))
}


if(colourBy=="foi"){
  coul <- colorRampPalette(base.pal)( length(unique(
                                                    METAPLOT[[1]][[colourBy]])
                                             ))
} else if(colourBy=="group"){
  #coul <- coul.group[match( unique(METAPLOT[[1]][[colourBy]]), 
  #                         group.v )
  #                   ]
  #coul <- coul.group[match( levels(METAPLOT[[1]][[colourBy]]), 
  #                          group.v )
  #                  ]
  coul.group <- c(colorRampPalette(base.pal)(length(featgrp.v)), "gray")
  coul.group <- adjustcolor(coul.group, alpha.f=0.5)
  names(coul.group) <- c(names(featgrp.v), "Rest") 
  coul <- coul.group[as.character(levels(METAPLOT$FCVSCP$group))]
}

foi.v.len <- length( unique(METAPLOT[[1]]$foi) )

if(foi.v.len < 15 | colourBy=="group"){
  num.col <- 1
  legend.tsize <- 20
} else {
  legend.tsize <- 5
  if(foi.v.len < 100){
    num.col <- 2
  } else {
    num.col <- 4
  }
}

plotParam <- list(geom_line(size=1, aes_string(colour=colourBy)
                            ),
                  geom_point(size=1.5, colour="black"),
                  scale_y_continuous(limits=c(-1,1)),
                  guides(colour=guide_legend(ncol=num.col)
                         ),
                  labs(y=bquote(bold("log"["2"]~"FC")), colour="Feature"),
                  bgr2, 
                  theme(legend.text=element_text(size=legend.tsize, 
                                                 face="bold"),
                        legend.title=element_text(size=12, face="bold"),
                        panel.grid.major.x=element_line(size=0.5, color="gray70",
                                                        linetype="dashed"))
                  )
#-------------------------------------------------------------------------------
# FCVSBIN plot
#-------------------------------------------------------------------------------
max.x <- max(METAPLOT$FCVSBIN$bin)
min.x <- min(METAPLOT$FCVSBIN$bin)

# Make legend label, foi with reference percbin
if(colourBy=="foi"){
  refpercbin <- paste0(round(METAPLOT$FCVSBIN[METAPLOT$FCVSBIN$bin==0, "percbin"], digits=2), "%")
  names(refpercbin) <- METAPLOT$FCVSBIN[METAPLOT$FCVSBIN$bin==0, "foi"]
  leg.lab <- paste(levels(METAPLOT[[1]][["foi"]]),
                   refpercbin[levels(METAPLOT[[1]][["foi"]])], sep="_")
}

p <- ggplot(data=METAPLOT$FCVSBIN, aes(x=bin, y=log2(percbinDIVref), group=foi)
            ) +
  scale_colour_manual(values=coul, labels=if(colourBy=="foi"){leg.lab}else{ levels(METAPLOT[[1]][["group"]]) }) +
  scale_x_continuous(labels=as.vector(rbind("",seq(from=min.x, to=max.x, by=2)))[-1],
                     breaks=min.x:max.x,
                     limits=c(min.x-0.5, max.x+0.5)
  ) +
  labs(title=paste0(id, "_order=pos", FCBIN.order), x="Contiguous 40-kb bins") + 
  plotParam

ggsave(filename=paste0(out.dir, "/", id, "_colBy", colourBy,
                       "_FCVSBIN_metaPcomb.pdf"), 
       width=15, height=15, plot=p)

# To keep colour of feature constant between plots
coul.ordrd <- match(levels(METAPLOT[[2]][[colourBy]]), 
                    levels(METAPLOT[[1]][[colourBy]])
                    )
#-------------------------------------------------------------------------------
# FCVSCP plot
#-------------------------------------------------------------------------------
# Make legend label, foi with reference percbin
if(colourBy=="foi"){
  refpercbin <- paste0(round(METAPLOT$FCVSCP[METAPLOT$FCVSCP$Cp==1, "percbin"], digits=2), "%")
  names(refpercbin) <- METAPLOT$FCVSCP[METAPLOT$FCVSCP$Cp==1, "foi"]
  leg.lab <- paste(levels(METAPLOT[[2]][["foi"]]),
                   refpercbin[levels(METAPLOT[[2]][["foi"]])], sep="_")
} 

p <- ggplot(data=METAPLOT$FCVSCP, aes(x=Cp, y=log2(percbinDIVref), group=foi)
            ) +
  scale_colour_manual(values=coul[coul.ordrd], labels=if(colourBy=="foi"){leg.lab}else{ levels(METAPLOT[[2]][["group"]]) }) + 
  scale_x_continuous(breaks=seq(from=1, to=21, by=1),
                     labels=as.vector(rbind(seq(from=1, to=21, by=2),
                                  rep(""))
                     )[-22]
                     ) +
  labs(title=paste0(id, "_order=Cp", FCCP.order), x=expression(bold( "c"["p"]) )
       ) + 
  plotParam 

ggsave(filename=paste0(out.dir, "/", id, "_colBy", colourBy,
                       "_FCVSCP_metaPcomb.pdf"), 
       width=15, height=15, plot=p)

# rm(list=ls()); gc()



