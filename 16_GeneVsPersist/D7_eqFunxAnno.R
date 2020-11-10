################################################################################
# Top 5% genes based on Cp/Cs and number of overlapping contacts for a given Cp/Cs
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start_time <- Sys.time()

whorunsit = "LiezelLinuxDesk" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
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
out.dir = paste0(objective.dir, "/out_/funxAnno")
### OTHER SETTINGS #############################################################
# 2(2MB gap) or "05"(0.5 MB minimum gap), refers to minimum gap accepted to classify a contact, 
# two points should be far enough to filter for contacts within a TAD
gcb = "min2Mb" #min05Mb
chr = "chrALL"
id.v = c("cp", "Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", 
         "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
nCPU = 2L
refseq = "ALL" 
anno.nme = "hg19anno"
# TRCOUNT.MX
src.affix = "ContactCountPerTr" 
Lref.v = c("orig", "ave", "m2sd")
# Regenerate data for GO functional analyses?  
regenerateData = TRUE  
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(ggplot2)
library(itertools)
source(paste0(lib, "/enGO.R"))
source(paste0(lib, "/enKEGG.R"))
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
annot <- fread( file=paste0(annofile.dir, "/", anno.nme, "_ALL"), 
                header=TRUE, data.table=FALSE, stringsAsFactors=FALSE,
                select=c("name2", "name", "txStart", "txEnd", "uniqueID"),
                colClasses=list(character=c("name2", "name", "uniqueID"),
                                integer=c("txStart", "txEnd")) )

annot$name <- gsub(x=annot$name, replacement="", pattern="[_].*$")
annot$txSize <- annot$txEnd-annot$txStart
# 1st reordering, by accession number(name), NM then NR
# 2nd reordering, decreasing txSize (hence the negative sign) so that first 
# transcript for a gene will be longest
annot <- annot[order(annot$name, -(annot$txSize), decreasing=FALSE), ]

id.v.len <- length(id.v)
toExport <- c("chr", "gcb", "refseq", "out.dir", "id.v", "regenerateData",
              "src.affix", "annot")
  
foreach( itr=isplitVector(1:id.v.len, chunks=nCPU), .inorder=FALSE, .combine="rbind",
         .export=toExport, .noexport=ls()[!ls()%in%toExport]
) %op% {
  
  for(i in itr){
    
    id <- id.v[i]
    affix <- paste0(chr, "_", gcb, "_", refseq, "_", id)
    
    if(regenerateData==TRUE){
      
      lst <- foreach(x=Lref.v, .inorder=TRUE, .combine="c") %do% {
        # Load TRCOUNT.MX
        load(file=paste0(objective.dir, "/out_/", affix, "_", src.affix, "_", 
                         x, ".RData"))
        
        # Remove irrelevant Cp/Cs values
        TRCOUNT.MX <- TRCOUNT.MX[, colSums(TRCOUNT.MX)!=0 & colnames(TRCOUNT.MX)!="0"]
        
        # Only consider transcripts overlapping with at least one contact
        TRCOUNT.MX <- TRCOUNT.MX[!rowSums(TRCOUNT.MX)==0,] 
        # Background genes
        uniqueID.v <- as.numeric(rownames(TRCOUNT.MX))
        # Get unique HUGO symbols of genes overlapping with contacts
        bgr <- unique(annot[uniqueID.v,"name2"])
        bgr.len <- length(bgr)
        # Top 5% relative to length of background genes
        genelist.len.ref <- floor(bgr.len*0.05)
        
        idVal <- as.character(
          sort( as.numeric(colnames(TRCOUNT.MX)), decreasing=TRUE )
        )
        idVal.len <- length(idVal)
        genelist <- NULL; genelist.len <- 0; cnt <- 1; out <- NULL
        
        while(genelist.len < genelist.len.ref){
          test <- TRCOUNT.MX[,idVal[cnt]]>0
          count.v <- TRCOUNT.MX[ test,idVal[cnt] ]
          uniqueID.v <- rownames(TRCOUNT.MX)[test]
          rm(test)
          df <- annot[annot$uniqueID%in%uniqueID.v,c("name2", "uniqueID", "txSize", "name")]
          ind <- match(df$uniqueID, uniqueID.v)
          df <- data.frame(df, idVal=rep(idVal[cnt]), count=count.v[ind], 
                           stringsAsFactors=FALSE)
          rm(count.v, uniqueID.v)
          
          addout <- by(data=df[, c("uniqueID", "txSize", "idVal", "count")], 
                       INDICES=df$name2, FUN=function(x){
                         
                         # Select transcript with maximum count of contact overlap
                         # If > 1 transcript, select the first entry in annotation file;
                         # should be the longest (or one of) NM transcript (or NR if no 
                         # NM among the longest) because annot was reordered
                         # prioritising NM
                         ind <- which(x$count==max(x$count))
                         if(length(ind)>1L){ind<-ind[1]}
                         x[ind,]
                       })
          addgenes <- setdiff(names(addout), genelist)
          genelist <- c(genelist, addgenes)
          genelist.len <- length(genelist)
          
          addout <- setDT(x=do.call( "rbind", addout[addgenes] ), 
                          keep.rownames="name2")
          if(nrow(addout)!=0){ addout <- addout[order(addout$count, decreasing=TRUE),] }
          out <- rbind(out, addout)
          
          rm(addout, addgenes, df);gc()
          
          print(paste0(idVal[cnt], ",", genelist.len))
          
          cnt <- cnt + 1L
        }
        
        out.len <- nrow(out)
        if(out.len > genelist.len.ref){
          # Trim out data table from the bottom, still satisfying genelist.len.ref
          # Rev because bottom groups should come first 
          group.id <- rev( paste(out$idVal, out$count, sep=",")  )
          agg <- table(group.id)
          # Reorder based on out data table
          agg <- agg[ order( match(names(agg), unique(group.id)) ) ]
          
          toSub <- out.len-genelist.len.ref
          agg.len <- length(agg)
          agg <- sapply(X=1:agg.len, FUN=function(x){
            sum(agg[1:x])
          })
          
          y <- agg[agg<=toSub]
          if( length(y)!=0 ){ out <- out[-( (out.len-max(y)+1L):out.len ),] }
          rm(agg, group.id, agg.len); gc()
        }
        
        print(x)
        
        return(list(out, bgr))
        
      } # foreach loop end
      
      # Name lst
      prefix <- paste0(">unique_genes_", id, "_")
      suffix.v <- c("_top5_end", "_bgr_end")
      nme <- list()
      for(l in Lref.v){
        nme[[l]] <- paste(prefix, l, suffix.v, sep="")
      }
      names(lst) <- unlist(x=nme, use.names=FALSE)
      
      # Save data on top 5 percent genes
      ind.genes <- grep(x=names(lst), pattern="top5", fixed=TRUE)
      T5PGENES <- lst[ind.genes]
      save(T5PGENES, file=paste0(out.dir, "/", affix, "_top5percGenes.RData"))
      rm(T5PGENES); gc()
      
      # Retrieve gene names for making text files for GO analyses
      top5genes <- sapply(X=Lref.v, simplify=FALSE, FUN=function(dt.nme){
        dt <- lst[[paste0(">unique_genes_", id, "_", dt.nme, "_top5_end")]]
        genes <- dt[,"name2"][[1]]
      })
      
      # Generate text file with the selected genes with and without length 
      # equalisation and the background to be used for GO analysis which 
      # is all genes overlapping with HiC contacts
      lst[ind.genes] <- top5genes; rm(top5genes)
      
      # Tally number of genes in genelist and background
      info <- paste( c( id, unlist(lapply(X=lst, FUN=length)) ), collapse=";" )
      # G=genelist, B=background
      write(x=info, file=paste0(out.dir, "/", chr, "_", gcb, "_", 
                                refseq, "_GB", paste(Lref.v, collapse="")),
            append=TRUE)
      
      GO.data <- unlist(
        
        rbind(
          names(lst),
          lapply(X=lst, FUN=function(x){ paste(x, collapse=";") })
        )
        
      )
      
      write(x=GO.data, file=paste0(out.dir, "/", chr, "_", gcb, "_", 
                                   refseq, "_top5genes"),
            append=TRUE)
      rm(lst, info); gc()
      
    } else {
      GO.data <- readLines(con=paste0(out.dir, "/", chr, "_", gcb, "_", 
                                      refseq, "_top5genes"))
      GO.data <- grep( x=GO.data, pattern=paste0(">unique_genes_", id, "_") )
    } 
    
    # HUGOsymbol-entrezID conversion table
    tbl <- fread(file=paste0(annofile.dir, "/", anno.nme, "_uniqueHUGOandEntrezID"),
                 data.table=FALSE, header=TRUE, stringsAsFactors=FALSE)  
    
    for( i in Lref.v ){
      
      genelist <- grep(x=GO.data, pattern=paste0(i, "_top5"), fixed=TRUE)
      genelist <- strsplit(x=GO.data[genelist+1L], split=";")[[1]]
      bgr <- grep(x=GO.data, pattern=paste0(i, "_bgr"), fixed=TRUE)
      bgr <- strsplit(x=GO.data[bgr+1L], split=";")[[1]]
      
      ego <- enGO(genes=genelist,
                  # Organism Db object
                  orgDb=org.Hs.eg.db,
                  # c("SYMBOL", "ncbi-geneid")
                  inputKey="SYMBOL",
                  ontology=c("BP","CC","MF"),
                  # Background genes
                  bgr=bgr,
                  # For saving table; if NULL, table not saved
                  filePath=paste0(out.dir, "/", affix, "_", i, "_GO"))
      
      gene.entrezIds <- tbl[match(x=genelist, table=tbl$hgnc_symbol), "entrezgene"]
      gene.entrezIds <- gene.entrezIds[!is.na(gene.entrezIds)]
      
      ekegg <- enKEGG(genes=gene.entrezIds,
                      org="hsa",
                      # c("kegg", "ncbi-geneid", "ncib-proteinid", "uniprot")
                      inputKey="ncbi-geneid",
                      # Background genes
                      bgr=bgr,
                      # For saving table; if NULL, table not saved
                      filePath=paste0(out.dir, "/", affix, "_", i, "_KEGG"))
    }
    
  } # itr for loop end
  
} # id.v for loop end

end_time <- Sys.time()
end_time-start_time  

# rm(list=ls())

  



