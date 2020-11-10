################################################################################
# Binary heatmap of GO/KEGG terms combining all tissues + cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start_time <- Sys.time()

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    objective.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    annofile.dir = "/Users/ltamon/Database/ucsc_tables/hsa_geneAnno"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    objective.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    annofile.dir = "/t1-data/user/ltamon/Database/ucsc_tables/hsa_geneAnno"
    os = "Linux"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    objective.dir = "/home/ltamon/DPhil/GenomicContactDynamics/4_AnnotationVsPersist"
    annofile.dir = "/home/ltamon/Database/ucsc_tables/hsa_geneAnno"
    os = "Linux"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
output.dir = paste0(objective.dir, "/out_eqFunxAnno")
### OTHER SETTINGS #############################################################
# 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
# two points should be far enough to filter for contacts within a TAD
gcb = "min2Mb" #min05Mb
chr = "chrALL"
refseq = "ALL" 
appro.v = c("GO_BP", "GO_MF", "GO_CC", "KEGG")
id.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
         "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC", "cp")
#id.v = c("cp", "Co", "Hi", "Lu")
Lref.v = c("orig", "ave", "m2sd")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(gplots)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
id.v.len <- length(id.v)
Lref.v.len <- length(Lref.v)

for(appro in appro.v){
  for(i in 1:id.v.len){
    
    id <- id.v[i]
    affix <- paste0(chr, "_", gcb, "_", refseq, "_", id)
    appro.lab <- paste0(appro, "_terms")
  
    lst <- lapply(X=Lref.v, FUN=function(x){
      
      df <- fread(file=paste0(output.dir, "/", affix, "_", x, "_", appro), 
            data.table=FALSE, header=TRUE, stringsAsFactors=FALSE,
            fill=TRUE, blank.lines.skip=FALSE, sep=",") 
      
      if( !all( c("p.adjust", "Description")%in%colnames(df) ) ){
        df <- data.frame( character(), integer() )
        colnames(df) <- c( appro.lab, paste0(id, "_", x) )
      } else {
        df <- df[,c("Description", "p.adjust")]
        colnames(df) <- c( appro.lab, paste0(id, "_", x) )
      }
      return(df)
    })
    names(lst) <- Lref.v
    tbl.temp <- lst[[1]]
    
    if(Lref.v.len!=1L){
      for(k in 2:Lref.v.len)
      tbl.temp <- merge(x=tbl.temp, y=lst[[k]], by=appro.lab, all=TRUE)
    } 
    
    if(i==1){
      mtbl <- tbl.temp
    } else {
      mtbl <- merge(x=mtbl, y=tbl.temp, by=appro.lab, all=TRUE)
    }
    
    print(x=paste0(appro, "_", id), quote=FALSE)
    
    rm(lst); gc()
    
  } # id.v.len for loop end
  
  # Make terms as rownames
  row.names(mtbl) <- mtbl[[1]]
  mtbl <- data.matrix(mtbl[,-1])
  
  # Make sure values are less than cutoff
  if( !all(mtbl[!is.na(mtbl)] < 0.05) ){
    stop("Not all values are less than 0.05 pvalue cutoff.")
  }
    
  mtbl[!is.na(mtbl)] <- 1L
  mtbl[is.na(mtbl)] <- 0L
  
  if(nrow(mtbl)==0){ next }
  
  pdf(file=paste0(output.dir, "/", chr, "_", gcb, "_", refseq,
                  "_", appro, "_funxannoBHM.pdf"), 
      width=10, height=10)
  heatmap.2(mtbl, Rowv=TRUE, Colv=FALSE, dendrogram="row", 
            scale="none", trace="none", na.rm=FALSE, na.color="darkblue", 
            margins=c(9, 7), col=c("gray80", "darkred"), cexRow=0.2, cexCol=0.5, 
            rowsep=seq(1:nrow(mtbl)), colsep=seq(from=Lref.v.len, to=ncol(mtbl), 
                                                 by=Lref.v.len), 
            sepwidth=c(0.01, 0.01), key=FALSE, 
            distfun = function(x) dist(x, method = "binary") )
  dev.off()
  
  rm(mtbl); gc()
  
} # appro.v for loop end

end_time <- Sys.time()
end_time-start_time  

# rm(list=ls())

  



