################################################################################
# Boxplot of minimum repeat count per Cp (combining all chromosomes)
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"
if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else if (whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/4_RepeatVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
rep.group = "subfam12" # "fam" | "subfam"
agerank.dir = paste0(wk.dir, "/Repeat_rankingbyAge")
minelm.dir = paste0(wk.dir, "/out_HicRepeatExploration/", rep.group)
out.dir = paste0(wk.dir, "/out_minRepCounts")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
# Number of repeat elements (372/56/62)
nCPU = 2L
# Age rank identifier
out.name = "GiorPubl_minRepCounts"
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
out.dir <- paste0(out.dir, "/", rep.group)
if( !dir.exists(out.dir) ){
  dir.create(out.dir)
}

print(paste0(gcb, " ", rep.group, "..."), quote=FALSE)

if(plotOnly==FALSE){
  
  col.nme <- ifelse(rep.group=="fam", "repFamily", "repName")
  agerank <- read.csv(file=paste0(agerank.dir, "/rep", rep.group, ".csv"),
                      header=TRUE, stringsAsFactors=FALSE)[,col.nme]; rm(col.nme)
  agerank.len <- length(agerank)
  
  ntis.v <- 1:21
  chr.v.len <- length(chr.v)
  for(ch in 1:chr.v.len){
    
    chr <- chr.v[ch]
    print(chr)
    load(paste0(minelm.dir,"/", chr, "_MinElm_", gcb, ".RData"))
    
    elements <- intersect(agerank, dimnames(MINELM.MX)[[2]][-1])
    elements.len <- length(elements)
    
    if(ch==1){
      
      # Initialize list
      lst <- vector("list", length(ntis.v))
      names(lst) <- as.character(ntis.v)
      MINREPCOUNTS <- rep(list(lst), agerank.len)
      names(MINREPCOUNTS) <- agerank
      
      toExport <- c("MINELM.MX", "elements")
      
      #### PARALLEL EXECUTION #########
      
      # List of lists
      lst.temp <- foreach( itr=isplitVector(1:elements.len, chunks=nCPU), 
                           .inorder=TRUE, .combine="c", .export=toExport, 
                           .noexport=ls()[!ls()%in%toExport]
      ) %op% {
        
        chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
          elm <- elements[i]
          print(elm)
          lst.ntis <- vector("list", 21)
          names(lst.ntis) <- 1:21
          
          for(ntis in 1:21){
            test <- MINELM.MX[,"ntis"]==ntis
            
            if(sum(test)!=0){
              lst.ntis[[ntis]] <- table(MINELM.MX[test,elm])
            } else {
              lst.ntis[[ntis]] <- NULL
            }
            
            print(ntis)
          } # ntis. v for loop end
          
          return(lst.ntis)
          
        }) # itr sapply end
        return(chunk)
      }
      
      ### END OF PARALLEL EXECUTION ###
      
      MINREPCOUNTS[elements] <- lst.temp
      rm(lst.temp); gc()
      
    } else {
      
      toExport <- c("MINELM.MX", "elements", "MINREPCOUNTS")
      
      #### PARALLEL EXECUTION #########
      
      # List of lists
      lst.temp <- foreach( itr=isplitVector(1:elements.len, chunks=nCPU), 
                           .inorder=TRUE, .combine="c", .export=toExport, 
                           .noexport=ls()[!ls()%in%toExport]
      ) %op% {
        
        chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
          elm <- elements[i]
          print(elm)
          
          lst.ntis <- vector("list", 21)
          names(lst.ntis) <- 1:21
          for(ntis in 1:21){
            test <- MINELM.MX[,"ntis"]==ntis
            
            if(sum(test)!=0){
              
              v <- c(MINREPCOUNTS[[elm]][[ntis]], table(MINELM.MX[test,elm]))
              agg <- aggregate(x=v, by=list(names(v)), FUN=sum)
              agg.v <- agg$x
              names(agg.v) <- agg$Group.1
              lst.ntis[[ntis]] <- agg.v
              
            } else {
              lst.ntis[[ntis]] <- NULL
            }
            print(ntis)
          } # ntis. v for loop end
          
          return(lst.ntis)
          
        }) # itr sapply end
       return(chunk)
      }
      
      ### END OF PARALLEL EXECUTION ###
      
      MINREPCOUNTS[elements] <- lst.temp
      rm(lst.temp); gc()
      
    }
    
    rm(MINELM.MX, elements); gc()
    
  } # chr.v for loop end
  
  save(MINREPCOUNTS, 
       file=paste0(out.dir, "/chrALL_", gcb, "_", out.name, ".RData"))
  
} else {
  
  # Load MINREPCOUNTS
  load(file=paste0(out.dir, "/chrALL_", gcb, "_", out.name, ".RData"))

}

drop <- unlist(lapply(MINREPCOUNTS, FUN=length))
MINREPCOUNTS <- MINREPCOUNTS[drop==21]
elements <- names(MINREPCOUNTS)

affix <- paste0("chrALL_", gcb, "_", out.name)

for(elm in elements){
  
  lst <- MINREPCOUNTS[[elm]]
  df <- sapply(X=names(lst), simplify=FALSE, 
               FUN=function(ntis){
                 v <- lst[[ntis]]
                 mincount <- as.numeric( names(lst[[ntis]]) )
                 names(v) <- NULL
                 # mapply(FUN, element to repeat, number of times)
                 df <- cbind(ntis=rep( as.numeric(ntis) ), 
                             mincount=unlist( mapply(rep, mincount, v) ) 
                             )
               })
                 
  df <- do.call("rbind", df)
  colnames(df) <- c("cp", "mincount")
  affix1 <- paste0( affix, "_", gsub(pattern="[^[:alnum:][:space:]]", 
                                    replacement="", x=elm) )
  
  png(file=paste0(out.dir, "/", affix1, "_bp.png"), 
      res=300, width=3000, height=3000)
  boxplot(formula=mincount~cp, data=df, col="honeydew3", 
          cex.lab=1.3, cex.axis=1.3,
          outline=FALSE, xlab="", ylab="", main="")
  # x axis
  mtext(side=1, text=expression("c"["p"]), line=3, cex=1.5)
  # y axis
  mtext(side=2, text="Contact min. repeat count", line=2.7, cex=1.5)
  # Diagram title
  mtext(side=3, text=affix1, line=1.5, cex=1.5)
  dev.off()

  rm(df)
  
  print(elm)
  
}
  
# rm(list=ls()); gc()