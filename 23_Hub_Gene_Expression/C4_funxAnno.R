################################################################################
# Functional enrichment analysis of genes overlapping with contacting regions 
# forming hubs
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/SahakyanLab/CoreGenomeExplorer"
    data.dir= "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
hub.dir = paste0(wk.dir, "/out_hubfile")
out.dir = paste0(wk.dir, "/out_funxAnno/bgrPerChr")
annofilePath = paste0(wk.dir, "/txTable/hg19anno_ALL")
# Converter of HUGO gene symbols to ncbi-geneid for KEGG
hugoEntrezPath = paste0(wk.dir, "/txTable/hg19anno_SYMBOLtoENTREZID_052020")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
hub.id = "FC_topCP3_gapBin50"
plotOnly = FALSE
bgrPerChr = TRUE # If FALSE, uses all unique HUGO symbols
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(org.Hs.eg.db)
ORG = c(GO="org.Hs.eg.db", KEGG="hsa")
library(clusterProfiler)
library(gplots)
source(paste0(lib, "/funxAnno.R"))
#source(paste0(lib, "/convertGeneKeyType.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(plotOnly==FALSE){
  
  bgrgenes <- read.delim(file=annofilePath, stringsAsFactors=FALSE, 
                         header=TRUE)[,c("chrom", "name2")]
  bgrgenes <- bgrgenes[!duplicated(bgrgenes$name2),]
  writeLines(bgrgenes[,"name2"], con=paste0(out.dir, "/hg19_bgrALL.txt"))
  
  HE.conv <- read.delim(file=hugoEntrezPath, header=TRUE, row.names=NULL, 
                        stringsAsFactors=FALSE)
  HE.conv <- HE.conv[!is.na(HE.conv$SYMBOL) & !is.na(HE.conv$ENTREZID),]
  #HE.conv <- convertGeneKeyType(genes=NULL, convTablePath=hugoEntrezPath, drop.NA=TRUE)
  HE.conv$ENTREZID <- as.character(HE.conv$ENTREZID)
  bgrgenes.KEGG <- HE.conv[HE.conv$SYMBOL%in%bgrgenes$name2,]
  bgrgenes.KEGG <- bgrgenes.KEGG[!duplicated(bgrgenes.KEGG$ENTREZID),]
  bgrgenes.KEGG <- merge(x=bgrgenes.KEGG, y=bgrgenes, by.x="SYMBOL",
                         by.y="name2", all.x=TRUE)
  writeLines(bgrgenes.KEGG[,"ENTREZID"], con=paste0(out.dir, "/hg19_bgrALL_entrezid.txt"))
  
  if(!bgrPerChr){
    chr.bgr.TF <- rep(TRUE, times=nrow(bgrgenes))
    chr.bgrKEGG.TF <- rep(TRUE, times=nrow(bgrgenes.KEGG))
    bgr.id <- "bgrAll"
  } 
  
  for(chr in chr.v){
    
    if(bgrPerChr){
      chr.bgr.TF <- bgrgenes$chrom==chr
      chr.bgrKEGG.TF <- bgrgenes.KEGG$chrom==chr
      bgr.id <- "bgrPerChr"
    } 
    
    out.id <- paste0(gcb, "_", chr, "_", hub.id)
    hub.v <- list.files(path=hub.dir, pattern=out.id, full.names=FALSE) 
    if(length(hub.v)==0){
      print(paste0("No hubs for ", out.id), quote=FALSE)
      rm(out.id, hub.v)
      next
    }
    hub.v.len <- length(hub.v)
    
    for(h in 1:hub.v.len){
      
      hub <- hub.v[h]
      fgrgenes <- read.csv(file=paste0(hub.dir, "/", hub), header=TRUE, 
                           stringsAsFactors=FALSE)[,"gene"]
      if( all(is.na(fgrgenes)) ){
        print(paste0("No gene for ", hub, "."), quote=FALSE)
        next
      }
      fgrgenes <- unlist(strsplit(x=fgrgenes, split=";"))
      fgrgenes <- unique(fgrgenes[!is.na(fgrgenes)])
      
      # Functional annotation
      funxAnnoOut <- sapply(X=c("GO_ALL", "KEGG"), simplify=FALSE, FUN=function(appro){
        if(appro=="KEGG"){
          # Convert fgrgenes to ENTREZID
          fgrgenes.KEGG <- unique(HE.conv[HE.conv$SYMBOL%in%fgrgenes, "ENTREZID"])
          bgr <- unique(bgrgenes.KEGG[chr.bgrKEGG.TF,"ENTREZID"])
          writeLines(bgr, con=paste0(out.dir, "/hg19_", chr, "_bgr_entrezid.txt"))
          df <- funxAnno(input=list(fgrgenes.KEGG, bgr), 
                         org=ORG["KEGG"], inputKey="ncbi-geneid", approach=appro,
                         filePath=NULL)
          rm(fgrgenes.KEGG, bgr)
          # Convert geneID in output to HUGO symbol
          df$geneID <- sapply(X=df$geneID, simplify=TRUE, FUN=function(id){
            id <- unique(strsplit(x=id, split="/", fixed=TRUE)[[1]])
            paste(x=HE.conv[HE.conv$ENTREZID%in%id,"SYMBOL"],
                  collapse="/")
          })
        } else {
          bgr <- unique(bgrgenes[chr.bgr.TF,"name2"])
          writeLines(bgr, con=paste0(out.dir, "/hg19_", chr, "_bgr.txt"))
          df <- funxAnno(input=list(fgrgenes, bgr),
                         org=ORG["GO"], inputKey="SYMBOL", approach=appro, filePath=NULL)
          rm(bgr)
        }
        return(df)
      }); rm(fgrgenes)
      funxAnnoOut <- do.call("rbind", funxAnnoOut)
      rownames(funxAnnoOut) <- NULL
      hub <- strsplit(hub, split=".csv", fixed=TRUE)[[1]]
      write.csv(funxAnnoOut, file=paste0(out.dir, "/", hub, "_", bgr.id, ".csv"),
                row.names=FALSE, quote=FALSE)
      
      #-------------------Collect p.adjust for heatmap
      if(nrow(funxAnnoOut)==0){
        print(paste0("No terms for ", hub, "."), quote=FALSE)
        next
      }
      funxAnnoOut$Description <- paste0(funxAnnoOut$ONTOLOGY, ":", funxAnnoOut$Description)
      funxAnnoOut <- funxAnnoOut[,c("Description", "p.adjust")]
      colnames(funxAnnoOut) <- c("Terms", hub)
      if( !exists("hmap") ){
        hmap <- funxAnnoOut
      } else {
        hmap <- merge(x=hmap, y=funxAnnoOut, by="Terms", all=TRUE)
      }
      rm(funxAnnoOut); gc()
      print(paste0(hub, " done!"), quote=FALSE)
      
    } # hub.v.len for loop end
    
  } # chr.v for loop end
  
  if(exists("hmap")){
    save(hmap, file=paste0(out.dir, "/", hub.id, "_", bgr.id, "_hmap.RData"))
  } else {
    stop("hmap does not exist.")
  }
  
} else {
  fle.nme <- paste0(out.dir, "/", hub.id, "_", bgr.id, "_hmap.RData")
  if( file.exists(fle.nme) ){
    load(file=fle.nme)
  } else {
    stop("hmap file does not exist.")
  }
}

#-------------------Binary heatmap
for( appro in c("GO_BP", "GO_CC", "GO_MF", "KEGG") ){
  
  incl.TF <- grepl(x=hmap$Terms, pattern=paste0(appro, ":"), fixed=TRUE)
  if(sum(incl.TF)==0){ next }
  mx <- data.matrix(hmap[incl.TF,-1])
  dimnames(mx)[[1]] <- hmap[incl.TF,"Terms"]
  rm(incl.TF)

  zero.TF <- mx>=0.05 | is.na(mx)
  mx[zero.TF] <- 0
  mx[!zero.TF] <- 1
  rm(zero.TF)
  
  pdf(file=paste0(out.dir, "/", hub.id, "_", appro, "_", bgr.id, "_hmap.pdf"),
      width=50, height=50)
  heatmap.2(x=mx, Rowv=TRUE, Colv=TRUE, dendrogram="both", scale="none", 
            trace="none", na.rm=FALSE, na.color="black", margins=c(10, 10), 
            col=c("gray80", "darkred"), cexRow=0.2, cexCol=0.5, rowsep=1:nrow(mx),
            key=FALSE, distfun=function(x) dist(x, method="binary"),
            main=paste0(hub.id, "_", appro, "_", bgr.id))
  dev.off()
  rm(mx); gc()
  print(paste0(appro, " heatmap done!"), quote=FALSE)
  
}

# rm(list=ls())



