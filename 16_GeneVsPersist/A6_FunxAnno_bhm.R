################################################################################
# Make heatmap combining the gene enrichment analyses results for all Cps
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
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(wk.dir, "/out_FunxAnno")
out.dir = paste0(wk.dir, "/out_FunxAnno_bhm")
### OTHER SETTINGS #############################################################
# gcb + refseq
prefix = "min2Mb_ALL"
tableID.v = c("cp", "Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
              "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
funxAnnoMethod.v = c("GO_BP", "GO_MF", "GO_CC", "KEGG")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(gplots)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(appro in funxAnnoMethod.v){
  for(id in tableID.v){
  
    #funxAnnotable <- fread(file=paste0(data.dir, "/", prefix, "_", id, "_", 
    #                                   appro, "_funxAnnotable.csv"), 
    #                       data.table=FALSE, header=TRUE, stringsAsFactors=FALSE,
    #                       fill=TRUE, blank.lines.skip=FALSE, sep=",")
    funxAnnotable <- read.csv(paste0(data.dir, "/", prefix, "_", id, "_", 
                                     appro, "_funxAnnotable.csv"), check.names=FALSE,
                              fill=TRUE, stringsAsFactors=FALSE)
    funxAnnotable.row <- nrow(funxAnnotable)
    funxAnnotable.col <- ncol(funxAnnotable)
    # In case of 0 or 1 row-tables, skip and record
    if(funxAnnotable.row<2 | !funxAnnotable.col%in%c(8, 23)){
      
      outstring <- paste0(appro, "_", id, "_", funxAnnotable.row, "row_", 
                          funxAnnotable.col, "col")
      print(x=outstring, quote=FALSE)
      write(x=outstring, file=paste0(out.dir, "/", prefix, "_skipped_binaryHeatmap"),
            append=TRUE)
      next
    }
    
    # Make terms as rownames
    row.names(funxAnnotable) <- funxAnnotable[[1]]
    funxAnnotable <- data.matrix(funxAnnotable[,-1])
    
    # Make sure values are less than cutoff
    if( !all(funxAnnotable[!is.na(funxAnnotable)] < 0.05) ){
      stop("Not all values are less than 0.05 pvalue cutoff.")
    }
    
    funxAnnotable[!is.na(funxAnnotable)] <- 1L
    funxAnnotable[is.na(funxAnnotable)] <- 0L
    
    pdf(file=paste0(out.dir, "/", prefix, "_", id, "_", 
                    appro, "_funxAnnotable.pdf"), width=10, height=7)
    heatmap.2(funxAnnotable, Rowv=TRUE, Colv=FALSE, dendrogram="row", 
              scale="none", trace="none", na.rm=FALSE, na.color="darkblue", 
              margins=c(9, 7), col=c("gray80", "darkred"), cexRow=0.2, cexCol=0.5, 
              rowsep=seq(1:nrow(funxAnnotable)), sepwidth=c(0.01, 0.01),
              key=FALSE, distfun = function(x) dist(x, method = "binary") )
    dev.off()
    
    print(x=paste0(appro, "_", id), quote=FALSE)
  }
  
}

# rm(list=ls())
