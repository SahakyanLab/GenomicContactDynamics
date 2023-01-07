################################################################################
# HicRepeatVsPersist heatmap
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"
if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = paste0("/Users/ltamon")
    wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")
  } else if (whorunsit == "LiezelCluster"){
    home.dir = paste0("/project/sahakyanlab/ltamon")
    wk.dir = paste0(home.dir, "/DPhil/GenomicContactDynamics/4_RepeatVsPersist")
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = paste0(home.dir, "/DPhil/lib")

rep.group = "subfam" # "fam" | "subfam" | "subfam6"
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
PreElmTissDyn.dir = paste0(wk.dir, "/out_HicRepeatHeatmapData/subfam_sumrep_atleast2sumrep") #, rep.group)
# hmclustPth = paste0(wk.dir, "/out_HicRepeatHeatmap/hm_famVssubfam_clust.csv")
out.dir = paste0(wk.dir, "/out_HicRepeatHeatmap/subfam_sumrep_atleast2sumrep")
### OTHER SETTINGS #############################################################
# Age rank identifier
out.name = "GiorPubl" #"GiorPubl" 
gcb = "min2Mb"
chr = "chrALL" 
# Regenerate ELMTISSDYN?
regenerateData = TRUE #FALSE
# Lineplot per repeat of average minimum repeat count per Cp
lineplot = TRUE
cluster.TF = FALSE
addLoess = TRUE # Only works if cluster.TF=FALSE
ylim.loess = c(-1,1)
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
library(viridis)
library(data.table)
library(ComplexHeatmap)
source(paste0(lib, "/doCorTest.R"))
source(paste0(wk.dir, "/lib/myheatmap.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * 
################################################################################
#out.dir <- paste0(out.dir, "/", rep.group)
#if( !dir.exists(out.dir) ){
#  dir.create(out.dir)
#}
linep.dir <- paste0(out.dir, "/lineplot")
cor.dir <- paste0(out.dir, "/correlation")
if( !dir.exists(linep.dir) & lineplot){ dir.create(linep.dir) }
if( !dir.exists(cor.dir) & addLoess ){ dir.create(cor.dir) }

#coul <- c(viridis::viridis(n=20)[c(4,7,9,18,20)], viridis::magma(n=20)[c(20:16)])
coul <- viridis(n=299)
#coul <- c(viridis(n=5)[1:4], 
#          # From viridis::plasma
#          c("#F0F921FF", "#FDB130FF", "#FCA338FF", "#E56A5DFF", "#DE6065FF"))

id <- paste0(chr, "_", gcb)
ELMTISSDYN <- list()
if(regenerateData){
  #---------------------------------------
  # Generate ELMTISSDYN.MX from PREELMTISSDYN.MX
  col.nme <- ifelse(rep.group=="fam", "repFamily", "repName")
  agerank <- read.csv(file=paste0(agerank.dir, "/rep", rep.group, ".csv"),
                      header=TRUE, stringsAsFactors=FALSE)[,col.nme]
  # Load PREELMTISSDYN.MX
  load(file=paste0(PreElmTissDyn.dir, "/", chr, "_", gcb, "_", out.name, 
                   "_PreElmTissDyn.RData"))
  num.contacts.each.ntis <- PREELMTISSDYN.MX[1,]
  ELMTISSDYN[["raw"]] <- t( apply(X=PREELMTISSDYN.MX[-1,], MARGIN=1, FUN=function(row){
    # 100*(row/num.contacts.each.ntis)
    row/num.contacts.each.ntis
  }) )
  
  ## Remove centr family
  #ELMTISSDYN$raw <- ELMTISSDYN$raw[ !rownames(ELMTISSDYN$raw) %in% c("centr"), ]
  # c("centr", "ERV", "scRNA", "srpRNA", "Other",
  #"SINE", "Helitron", "TcMar", "rRNA", "snRNA")
  #---------------------------------------
  # # Lineplot per repeat of average minimum repeat count per Cp
  # ntis <- 1:21
  # if(lineplot){
  #   element <- rownames(ELMTISSDYN$raw)
  #   element.len <- length(element)
  #   for(elm.ind in 1:element.len){
  #     elm <- element[elm.ind]
  #     vals <- ELMTISSDYN$raw[elm,]
  #     
  #     affix1 <- gsub(pattern="[^[:alnum:][:space:]]", replacement="", x=elm)
  #     affix1 <- paste0(affix1, "_", elm.ind)
  #     
  #     pdf(file=paste0(out.dir, "/lineplot/", id, "_lineplot_", affix1, ".pdf"), 
  #         width=10, height=10)
  #     plot(x=ntis, y=vals, col="#55bde6", cex.lab=1.3, cex.axis=1.3, 
  #          pch=19, lwd=5, xlab="", ylab="", main="")
  #     lines(x=ntis, y=vals, col="#55bde6")
  #     # X axis
  #     mtext(side=1, text=expression("c"["p"]), line=3, cex=1.5)
  #     # Y axis
  #     mtext(side=2, text="Contact fraction with non-0 min repeat count", 
  #           line=2.7, cex=1)
  #     # Diagram title
  #     mtext(side=3, text=paste0(id, "_", elm), line=1.5, cex=1.5)
  #     dev.off()
  #   } # element for loop end
  # }
  # pdf(file=paste0(linep.dir,  "/", id, "_densplot_all.pdf"), 
  #     width=10, height=10)
  # d <- density(ELMTISSDYN$raw, na.rm=T)
  # plot(d, xlab="", ylab="", main="", col="#55bde6")
  # polygon(d, col="#55bde6", border="#55bde6")
  # # X axis
  # mtext(side=1, text="values",#bquote(bold("log"["10"]~"metric")),
  #       line=3, cex=1.5)
  # # Y axis
  # mtext(side=2, text=expression(bold("Density")), line=2.7, cex=1.5)
  # # Diagram title
  # mtext(side=3, text=paste0(id, "_all_ELMTISSDYN$raw_N=", length(rownames(ELMTISSDYN$raw))),
  #       line=1.5, cex=1.5)
  # dev.off()
  #---------------------------------------
  # Normalize (z-score) values in ELMTISSDYN$raw for making heatmap
  ELMTISSDYN[["norm"]] <- ELMTISSDYN$raw
  for(i in 1:dim(ELMTISSDYN$raw)[1]){
    ELMTISSDYN$norm[i,] <- (ELMTISSDYN$raw[i,] - mean(ELMTISSDYN$raw[i,]))/
      sd(ELMTISSDYN$raw[i,])
  }
  ELMTISSDYN[["fc"]] <- log2(ELMTISSDYN$raw/ELMTISSDYN$raw[,"1"])
  save(ELMTISSDYN, file=paste0(out.dir, "/", chr, "_", gcb, "_ElmTissDyn_",
                               out.name, ".RData"))
} else {
  load(file=paste0(out.dir, "/", chr, "_", gcb, "_ElmTissDyn_", out.name, 
                   ".RData"))
}

# Generate heatmap
mx.nme.v <- names(ELMTISSDYN)
for(mx.nme in mx.nme.v){
  mx <- ELMTISSDYN[[mx.nme]]
 
  drp.rw <- apply( X=mx, MARGIN=1, FUN=function(rw) any(!is.finite(rw)) )
  if(mx.nme=="fc" & rep.group=="fam"){
    drp.rw[names(drp.rw)=="centr"] <- TRUE
  }
  mx[!is.finite(mx)] <- NA
  
  if(mx.nme=="raw"){
    at.v <- seq(from=0, to=1, by=0.1)
  } else if( mx.nme%in%c("norm", "fc") ){
    
    #if(mx.nme=="fc" & rep.group=="fam"){ 
    #  at.v <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 6, 8)
    #}  else {
      at.v <- seq(from=floor(min(as.vector(mx[!drp.rw]), na.rm=TRUE)), 
                  to=ceiling(max(as.vector(mx[!drp.rw]), na.rm=TRUE)), 
                  length.out=10)
      at.v <- sort(unique(c(0, round(x=at.v, digits=1))))
    #}
      
  }
  #-------------------Generate heatmap
  pdf(file=paste0(out.dir, "/", id, "_", mx.nme, "_heatmap.pdf"), 
      width=10, height=10)
  
  if(cluster.TF){
    # Rows with non-finite values; exclude in clustering then append later to
    # clustered heatmap.
    h <- myheatmap(mx=mx[!drp.rw,], colScheme=coul, rep.group=rep.group, 
                   mx.nme=mx.nme, cluster=cluster.TF, at.v=at.v)
    
    if( sum(drp.rw)>0 ){
      
      drp.mx <- mx[drp.rw,]
      
      if( any(!is.na(drp.mx)) ){
        
        atdrp.v <- seq(from=floor(min(as.vector(drp.mx), na.rm=TRUE)), 
                       to=ceiling(max(as.vector(drp.mx), na.rm=TRUE)), 
                       length.out=10)
        atdrp.v <- sort(unique(c(0, round(x=atdrp.v, digits=1))))
        
        
        
        if(mx.nme=="fc" & rep.group=="fam"){
          drp.mx <- rbind(mx["centr",], mx[drp.rw,])
          rownames(drp.mx) <- c("centr", rownames(mx[drp.rw,]))
          drp.mx <- drp.mx[!duplicated(rownames(drp.mx)),]
        }
        
        h2 <- myheatmap(mx=drp.mx, colScheme=coul, rep.group=rep.group, 
                        mx.nme=mx.nme, cluster=FALSE, at.v=atdrp.v)
        h <- h %v% h2
        
      } else {
        write(rownames(drp.mx), file=paste0(out.dir, "/", chr, "_", gcb, "_", out.name, 
                                            "_", mx.nme, "AllNArows_notinHmap.txt"))
        print(paste0(mx.nme, ": Rows with all NAs not included."), quote=FALSE)
      }
      
    }
    
  } else {
    h <- myheatmap(mx=mx, colScheme=coul, rep.group=rep.group, 
                   mx.nme=mx.nme, cluster=cluster.TF, at.v=at.v)
  }
  
  print(h)
  
  dev.off()
  
  rm(drp.rw)
  #-------------------Version of heatmap with loess
  if(addLoess){
    x.v <- as.numeric(colnames(mx))
    x.v.len <- length(x.v)
  
    # Loess line using weighted average Cp
    # loe <- apply(X=mx, MARGIN=1, FUN=function(rw){sum(x.v*rw)/x.v.len})

    # Loess line with correlation coefficient 
    
    element <- rownames(mx)
    element.len <- length(element)
    loe <- sapply(X=1:element.len, simplify=F, FUN=function(elm.ind){
      
      elm <- element[[elm.ind]]
      rw <- mx[elm.ind,]
      if( any(!is.finite(rw)) ){
        return(NaN)
      } else {
        
        cor.out <- doCorTest(xval=x.v, yval=rw, alt="two.sided", exactpval=F, 
                             out.dir=cor.dir, out.name=paste0(id, "_", mx.nme, "_", elm, "_cortest.RData"))
        
        cor.coef <- ifelse(cor.out$kend$p.value > 0.05, NA, unname(cor.out$kend$estimate))
        
      }
      
    })
    
    loe <- unlist(loe)
    names(loe) <- rownames(mx)
    
    info.log <- c(`Total element` = element.len, 
                  `With correlation data`=sum(is.finite(loe)),
                  `Not significant alpha 0.05`=sum(is.na(loe) & !is.nan(loe)),
                  `Missing value so not correlated`=sum(is.nan(loe)))
    info.log <- stack(info.log)
    colnames(info.log) <- c("count", "what")
    write.table(info.log, file=paste0(out.dir, "/", id, "_", mx.nme, "_loess_heatmap.txt"),
                row.names=F, quote=F, sep="\t")
    
    row.loe <- rowAnnotation(
      # `Weighted (z-score) \n average Cp`=
      #paste0("Spearman, >0.05 pval \n set to 0\n", perc.0.loe, "% 0val")=anno_lines(loe, smooth=TRUE, ylim=c(-1,1),
                                                                                  # gp=gpar(col="gray50", lwd=3), 
                                                                                  # pt_gp=gpar(col="black", cex=1.5),
                                                                                  # width=unit(2, "cm")),
      `Kendall, > 0.05 pval \n set to NA (no point)`=anno_lines(loe, smooth=TRUE, ylim=ylim.loess,
                                                                gp=gpar(col="gray50", lwd=3), 
                                                                pt_gp=gpar(col="black", cex=1.5),
                                                                width=unit(2, "cm")),
          
      annotation_name_gp=gpar(cex=0.5)
    )
    
    #h <- myheatmap(mx=mx[loe!=0,], colScheme=coul, rep.group=rep.group, 
    #               mx.nme=mx.nme, cluster=cluster.TF, at.v=at.v)
    
    pdf(file=paste0(out.dir, "/", id, "_", mx.nme, "_loess_heatmap.pdf"), 
        width=10, height=10)
    print(h+row.loe)
    dev.off()
    
    rm(row.loe)
  } # addLoess for loop end
  
  #---------------------------------------
  # Lineplot per repeat of average minimum repeat count per Cp
  ntis <- 1:21
  if(lineplot){
    element <- rownames(mx)
    element.len <- length(element)
    for(elm.ind in 1:element.len){
      elm <- element[elm.ind]
      vals <- mx[elm.ind,]
      
      if( !all(!is.finite(vals)) ){
        affix1 <- gsub(pattern="[^[:alnum:][:space:]]", replacement="", x=elm)
        affix1 <- paste0(affix1, "_", elm.ind)
        
        pdf(file=paste0(out.dir, "/lineplot/", id, "_lineplot_", affix1, "_", mx.nme, ".pdf"),
            width=10, height=10)
        plot(x=ntis, y=vals, col="#55bde6", cex.lab=1.3, cex.axis=1.3,
             pch=19, lwd=5, xlab="", ylab="", main="")
        lines(x=ntis, y=vals, col="#55bde6")
        # X axis
        mtext(side=1, text=expression("c"["p"]), line=3, cex=1.5)
        # Y axis
        mtext(side=2, text="Value in ELMTISSDYN",
              line=2.7, cex=1)
        # Diagram title
        mtext(side=3, text=paste0(id, "_", elm), line=1.5, cex=1.5)
        dev.off()
      }
  
    } # element for loop end
  }
  
  df <- as.data.frame(matrix(data=as.numeric(mx), ncol=1, dimnames=list(NULL, "value")))
  p <- ggplot(data=df, aes(x=value)) +
    geom_density(aes(y= ..scaled..), col="#55bde6", fill="#55bde6") +
    labs(y=NULL, title=paste0(id, "_allvaluesinELMTISSDYN_N=", length(mx[,1]))) + 
    bgr2
    
  ggsave(filename=paste0(linep.dir,  "/", id, "_densplot_all_", mx.nme, ".pdf"),
         width=10, height=10, units="in", plot=p)
  #---------------------------------------
  
  print(paste0(mx.nme, " done!"), quote=F)
  
  rm(mx, p, df)
  
}

# rm(list=ls()); gc()


