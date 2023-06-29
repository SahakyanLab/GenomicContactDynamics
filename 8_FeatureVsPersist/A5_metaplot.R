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

# Avoid left to right partial matching by $
options(warnPartialMatchDollar=TRUE)

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/8_FeatureVsPersist"
    data.dir = "/Users/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# Chromatin features directory
run.id = "inc"
foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/raw_associated")
foifile = paste0(wk.dir, "/foifile/FINAL_", run.id, "_foifile_priority_nperm10000_seed662_mxmskfr0_Cp21_pvalcutoff0.05_numOlapA_comOlap_edited_consistentAcrossCT")
#foifile = paste0(wk.dir, "/foifile/foifile_test")
# File of feature grouping
featgrpfile = paste0(wk.dir, "/features_group")
fetacp.dir = paste0(wk.dir, "/out_FETACP")
out.dir = paste0(wk.dir, "/out_metaplot")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
description = NULL
out.id = paste0("final_", run.id, "_priority_Cp21_numOlapAANDcomOlap_acrossTissues")
plotOnly = TRUE
# What bin to plot for FC vs. bin position
cp = 21
# What value to use for ordering?
FCBIN.order = 1 # pos 
FCCP.order = 21 # Cp

colourBy = "group" # "foi" | "group"
ylimits.v = c(-0.1, 1) # Watch out for missing points depending on limits
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(reshape)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
library(stringr)
library(yarrr); base.pal <- yarrr::piratepal("basel")
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/finaliseFOI.R"))
source(paste0(wk.dir, "/lib/finaliseFOI.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# File with feature groupings
featgrp <- readLines(con=featgrpfile)
header.ind <- grep(x=featgrp, pattern=">")
# Remove ">" from the headers
featgrp[header.ind] <- gsub(x=featgrp[header.ind], pattern=">", replacement="")
header.ind.len <- length(header.ind)
# List of patterns per group
featgrp.v <- sapply(X=1:header.ind.len, simplify=TRUE, FUN=function(x){
  
  first.pat <- header.ind[x]+1
  # NA if out of bounds
  last.pat <- ifelse(x==header.ind.len, length(featgrp), header.ind[x+1]-1)
  return( featgrp[first.pat:last.pat] )
  #return( paste(v, collapse="|") )

})

names(featgrp.v) <- featgrp[header.ind]
rm(featgrp)
feat.df <- stack(featgrp.v)
feat.df$ind <- as.character(feat.df$ind)

# Remove whitespaces
for( x in colnames(feat.df) ){
  feat.df[[x]] <- stringr::str_trim(string=feat.df[[x]], side="both")
  feat.df[[x]] <- stringr::str_squish(string=feat.df[[x]])
}

#-------------------------------------------------------------------------------

id <- paste0("chrALL_", gcb, "_desc_", out.id)

if(plotOnly==FALSE){
  
  # List of features
  foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)
  foi.v <- stringr::str_trim(string=foi.v, side="both")
  foi.v <- stringr::str_squish(string=foi.v)
  
  # Filter feature list by description
  if( !is.null(description) ){
    foi.v <- foi.v[ grepl(x=foi.v, pattern=paste0("_desc_", description), ignore.case=FALSE) ]
  }
  
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
    
    # Load FETACP.MX
    fetacp.id <- gsub(x=foi.v[f], pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
    load(file=paste0(fetacp.dir, "/chrALL_", gcb, "_", fetacp.id, "_fetacp.RData"))
  
    # Identify group of feature
    foipat <- paste0("foi_", strsplit(x=foi.v[f], split="foi\\_|\\_desc")[[1]][2], "_desc")
    grp <- feat.df$ind[ grep(x=paste0("foi_", feat.df$values, "_desc"), pattern=foipat) ]
    
    if( length(grp)>1 ){
      stop( paste0(fetacp.id, ": More than 1 group."))
    } else if( length(grp)==0 ){
      
      warning( paste0(fetacp.id, ": Can't be classified into a group; assigned to Unclassified"))
      grp <- "Unclassified"
      
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
    #df <- cbind(foi=rep(foi.v[f]), df)
    df <- cbind(foi=rep(fetacp.id), df)
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
    #df <- cbind(foi=rep(foi.v[f]), df)
    df <- cbind(foi=rep(fetacp.id), df)
    # Value for ordering plots
    ORDER$FCVSCP[f] <- log2(df[df$Cp==FCCP.order, "percbinDIVref"])
    FCVSCP[[ foi.v[f] ]] <- df
  
    # Add group
    FCVSBIN[[ foi.v[f] ]] <- cbind(group=rep(grp), FCVSBIN[[ foi.v[f] ]])
    FCVSCP[[ foi.v[f] ]] <- cbind(group=rep(grp),  FCVSCP[[ foi.v[f] ]])
    
    print(paste0(foi.v[f], " data obtained!"), quote=FALSE)
    
    rm(FETACP.MX, df, df.orig, meanCpref, meanCp1, fetacp.id, foipat)
    gc()
    
  } # foi.v for loop end
  
  METAPLOT <- list(FCVSBIN=do.call("rbind", FCVSBIN), 
                   FCVSCP=do.call("rbind", FCVSCP)) 
  # Order foi levels
  for( x in names(METAPLOT) ){
    
    rownames(METAPLOT[[x]]) <- NULL
    foi.ordrd <- order(ORDER[[x]], decreasing=TRUE)
    METAPLOT[[x]][["foi"]] <- factor(x=METAPLOT[[x]][["foi"]],
                                     levels=unique(METAPLOT[[x]][["foi"]])[foi.ordrd]
                                     )
    METAPLOT[[x]][["group"]] <- factor(x=METAPLOT[[x]][["group"]],
                                       levels=sort(unique(METAPLOT[[x]][["group"]]))
                                     )
    
  }
  
  save(METAPLOT, file=paste0(out.dir, "/", id, "_metaPcomb.RData"))
  
  rm(foi.ordrd); gc()
  
} else {
  # Load METAPLOT list
  load(file=paste0(out.dir, "/", id, "_metaPcomb.RData"))
}

if(colourBy=="foi"){
  
  foi.len <- length(levels(METAPLOT$FCVSBIN$foi))
  coul <- colorRampPalette(base.pal)(foi.len)
  
} else if(colourBy=="group"){
  
  group.len <- length(levels(METAPLOT$FCVSCP$group))
  coul.group <- c(colorRampPalette(base.pal)(group.len), "gray")
  coul.group <- adjustcolor(col=coul.group, alpha.f=0.7)
  names(coul.group) <- c(levels(METAPLOT$FCVSCP$group), "Unclassified") 
  coul <- coul.group[as.character(levels(METAPLOT$FCVSCP$group))]
  
}

foi.v.len <- length( unique(METAPLOT$FCVSBIN$foi) )

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
                  geom_hline(yintercept=0, linetype="dashed", colour="gray70", 
                             size=0.5),
                  scale_y_continuous(limits=ylimits.v
                                     ,breaks=seq(ylimits.v[1],ylimits.v[2],1)
                                     #,breaks=seq(0,1)
                                     ),
                  guides(colour=guide_legend(ncol=num.col)
                         ),
                  labs(y=bquote(bold("log"["2"]~"FC")), colour="Feature"),
                  bgr2, 
                  theme(legend.text=element_text(size=legend.tsize, 
                                                 face="bold"),
                        legend.title=element_text(size=12, face="bold"),
                        panel.grid.major.x=element_line(size=0.5, colour="gray70",
                                                        linetype="dashed"))
                  )
#-------------------------------------------------------------------------------
# FCVSBIN plot
#-------------------------------------------------------------------------------
max.x <- max(METAPLOT$FCVSBIN$bin)
min.x <- min(METAPLOT$FCVSBIN$bin)

# Make legend label, foi with reference percbin
if(colourBy=="foi"){
  
  refpercbin <- paste0( round(METAPLOT$FCVSBIN[METAPLOT$FCVSBIN$bin==0, "percbin"], digits=2), "%" )
  names(refpercbin) <- METAPLOT$FCVSBIN[METAPLOT$FCVSBIN$bin==0, "foi"]
  leg.lab <- paste(levels(METAPLOT$FCVSBIN[["foi"]]),
                   refpercbin[levels(METAPLOT$FCVSBIN[["foi"]])], sep="_")
  
}

p <- ggplot( data=METAPLOT$FCVSBIN, aes(x=bin, y=log2(percbinDIVref), group=foi) ) +
  scale_colour_manual(values=coul, 
                      labels=if(colourBy=="foi"){leg.lab}else{ levels(METAPLOT$FCVSBIN[["group"]]) }) +
  scale_x_continuous(labels=as.vector(rbind("",seq(from=min.x, to=max.x, by=2)))[-1],
                     breaks=min.x:max.x,
                     limits=c(min.x-0.5, max.x+0.5)
  ) +
  labs(title=paste0(id, "_order=pos", FCBIN.order), x="Contiguous 40-kb bins") + 
  plotParam

ggsave(filename=paste0(out.dir, "/", id, "_colBy", colourBy, "_FCVSBIN_metaPcomb.pdf"), 
       width=15, height=15, plot=p)

# To keep colour of feature constant between FCVSBIN and FCVSCP plots
coul.ordrd <- match(x=levels(METAPLOT$FCVSCP[[colourBy]]), 
                    table=levels(METAPLOT$FCVSBIN[[colourBy]])
                    )
#-------------------------------------------------------------------------------
# FCVSCP plot
#-------------------------------------------------------------------------------
# Make legend label, foi with reference percbin
if(colourBy=="foi"){
  
  refpercbin <- paste0(round(METAPLOT$FCVSCP[METAPLOT$FCVSCP$Cp==1, "percbin"], digits=2), "%")
  names(refpercbin) <- METAPLOT$FCVSCP[METAPLOT$FCVSCP$Cp==1, "foi"]
  leg.lab <- paste(levels(METAPLOT$FCVSCP[["foi"]]),
                   refpercbin[levels(METAPLOT$FCVSCP[["foi"]])], sep="_")
  
} 

p <- ggplot( data=METAPLOT$FCVSCP, aes(x=Cp, y=log2(percbinDIVref), group=foi) ) +
  scale_colour_manual(values=coul[coul.ordrd], 
                      labels=if(colourBy=="foi"){leg.lab}else{ levels(METAPLOT$FCVSCP[["group"]]) }) + 
  scale_x_continuous(breaks=seq(from=1, to=21, by=1),
                     labels=as.vector(rbind(seq(from=1, to=21, by=2),
                                  rep(""))
                     )[-22]
                     ) +
  labs( title=paste0(id, "_order=Cp", FCCP.order), x=expression(bold( "c"["p"])) ) + 
  plotParam 

ggsave(filename=paste0(out.dir, "/", id, "_colBy", colourBy, "_FCVSCP_metaPcomb.pdf"), 
       width=15, height=15, plot=p)

# rm(list=ls()); gc()



