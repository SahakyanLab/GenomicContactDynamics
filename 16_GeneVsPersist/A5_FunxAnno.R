################################################################################
# GO term and KEGG pathway enrichment for genes per Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
genelist.dir = paste0(wk.dir, "/out_anno_union")
out.dir = paste0(wk.dir, "/out_FunxAnno")
### OTHER SETTINGS #############################################################
# gcb + refseq
prefix = "min2Mb_ALL"
funxAnnoMethod.v = c("GO_MF", "GO_CC", "KEGG")
nCPU <- length(funxAnnoMethod.v)
contactFeature = c("cp", "count")
cutoffStr = 5L
celltiss.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
               "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
regenerateFunxAnnoTable = TRUE
# FALSE if you only want to save raw tables from functional annotation 
makePvalTable = TRUE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(org.Hs.eg.db)
library(clusterProfiler)
source(paste0(lib, "/UTL_doPar.R"))
source(paste0(lib, "/funxAnno.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
print(prefix, quote=FALSE)

# Generate forebackground strings defining foreground and background
foreback.combi <- c()
# Identifier of gene set
id.combi <- c()
if( !is.null(contactFeature) ){
  
  if("cp"%in%contactFeature){
    id.combi <- "cp"
    foreback.combi <- c( "_cp_1;all_genes", "_cp_1;_cp_HiC_all",
                         paste("_cp", paste(2:21, "_cp_1", sep=";"), sep="_") 
    )
  }
  
  if("count"%in%contactFeature){
    celltiss.v.len <- length(celltiss.v)
    id.combi <- c(id.combi, celltiss.v)
    # Basic unit
    mx <- rbind(
      c("count_1", "count_HiC_all"),
      cbind( paste("count", c(2:cutoffStr, paste0("grthan", cutoffStr)), sep="_"),
             rep("count_1"))
    )
    # Iterate over celltiss.v
    mx <- matrix(paste("", celltiss.v, rep(mx, each=celltiss.v.len), sep="_"), 
                 ncol=2)
    
    foreback.combi <- c(
      # add "<tissue>_count_1;all_genes_end"
      paste0("_", celltiss.v, "_count_1;all_genes"), 
      paste(mx[,1], mx[,2], sep=";"),
      foreback.combi
    )
    
    rm("mx", "celltiss.v", "celltiss.v.len"); gc()
  }
  
} else {
  stop("Invalid input for contactFeature (cp and/or count).")
}
#---------------------------------------
toExport <- c("genelist.dir", "out.dir", "prefix", "foreback.combi", "appro",
              "id.combi", "regenerateFunxAnnoTable", "makePvalTable")

#### PARALLEL EXECUTION #########

# Perform functional annotation
#for(appro in funxAnnoMethod.v){
foreach(appro=funxAnnoMethod.v, .inorder=FALSE, .export=toExport, 
        .noexport=ls()[!ls()%in%toExport]
) %op% {
  
  if(appro=="KEGG"){
    genes.str <- readLines(con=paste0(genelist.dir, "/", prefix, "_entrezID"))
    key <- "ncbi-geneid"
    org <- "hsa"
  } else if( appro%in%c("GO_BP", "GO_MF", "GO_CC") ){
    genes.str <- readLines(con=paste0(genelist.dir, "/", prefix, "_name2"))
    key <- "SYMBOL"
    org <- "org.Hs.eg.db"
  } else {
    stop("Invalid funxAnnoMethod input.")
  }
  lab <- paste0(appro, "_terms")
  
  for(id in id.combi){
    # cp here is treated as identifier along with celltiss.v, 
    # No "cp" in 21 tissues/cellline so no conflicts
    # Subset according to id
    
    # Id should be "_<id>_" to separate LC from TLC
    id <- paste0("_", id, "_")
    foreback.combi.sub <- foreback.combi[ grep(x=foreback.combi, 
                                               pattern=id, fixed=TRUE) ]
    foreback.combi.sub.len <- length(foreback.combi.sub)
    
    for(i in 1:foreback.combi.sub.len){
      lab1 <- gsub(x=foreback.combi.sub[i], pattern="\\;|\\_", 
                   replacement="")
      
      genes.lst <- sapply(X=strsplit(x=foreback.combi.sub[i], split=";")[[1]], simplify=FALSE, 
                          FUN=function(x){
                            x <- paste0(x, "_end")
                            ind <- grep(pattern=x, x=genes.str, fixed=TRUE)
                            if(length(ind)==1){
                              genes <- genes.str[ind+1L]; rm(ind)
                            } else {
                              stop("Pattern not selective.")
                            }
                            # Also split with comma because entrezIds (when gene names are converted to
                            # entrezIds for KEGG) of one gene are separated by comma
                            genes <- strsplit(x=genes, 
                                              split=paste0("\\;|\\,"))[[1]]
                            unique(genes[!is.na(genes) & genes!="NA" & genes!="" & genes!=" "])
                          })
      
      if(regenerateFunxAnnoTable==TRUE){
        funxAnnoOut <- funxAnno(input=genes.lst, org=org, inputKey=key, approach=appro, 
                                filePath=paste0(out.dir, "/", prefix, id, lab1, "_", 
                                                appro, ".csv"))[,c("Description", "p.adjust")]
      } else {
        funxAnnoOut <- read.csv(file=paste0(out.dir, "/", prefix, id, lab1, "_", appro, ".csv"),  
                                header=TRUE, stringsAsFactors=FALSE)[,c("Description", "p.adjust")]
      }
      rm(genes.lst)
      
      if(makePvalTable==TRUE){
        colnames(funxAnnoOut) <- c(lab, foreback.combi.sub[i])
        if(i==1){
          maintable <- funxAnnoOut
        } else {
          maintable <- merge(x=maintable, y=funxAnnoOut, by=lab, all=TRUE)
        }
      }
      rm(funxAnnoOut, lab1); gc()
      
    } # foreback.combi.sub.len for loop end
    
    if(makePvalTable==TRUE){
      write.csv(x=maintable, file=paste0(out.dir, "/", prefix, id, appro,
                                         "_funxAnnotable.csv"),
                row.names=FALSE, quote=FALSE)
      rm(maintable)
    }
    
    rm(foreback.combi.sub, foreback.combi.sub.len, id); gc()
    
  } # id.combi end for loop
  
} 

### END OF PARALLEL EXECUTION ###

# rm(list=ls())