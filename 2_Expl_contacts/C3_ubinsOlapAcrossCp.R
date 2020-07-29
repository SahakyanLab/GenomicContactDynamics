################################################################################
# Determine fraction overlap of unique contacting bins (long-range contact bins) 
# between Cp category. 
# query_ref, means that fraction will be calculated with respect to count of ref
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/2_Expl_contacts"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/2_Expl_contacts"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
out.dir = paste0(wk.dir, "/out_ubinsOlapAcrossCp")
### OTHER SETTINGS #############################################################
gcb = "min05Mb" #c("min2Mb", "min05Mb")
chr.v = paste0("chr", c(1:22, "X"))
Ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
         "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
Cp.v = 1:21
nCPU = 1L

query.v = c(as.character(1:21), "21", "1")
ref.v = c(rep(x="HiCAll", times=21), "1", "21")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(foreach)
library(doParallel)
library(itertools)
source(paste0(lib, "/UTL_doPar.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
ref.v <- as.character(ref.v)
query.v <- as.character(query.v)
query.v.len <- length(query.v)
if( query.v.len!=length(ref.v) ){ stop("Query and ref lengths not equal.") }

Cp.v <- sort(Cp.v, decreasing=FALSE)
ct.v <- c(Ct.v, "allCT"); ct.v.len <- length(ct.v)
chr.v.len <- length(chr.v)

chrLen.df <- read.delim(file=chrLenfile, stringsAsFactors=FALSE, header=TRUE)
totbin <- sum(chrLen.df[chrLen.df$chromosome%in%chr.v,"bins.40kb"])
rm(chrLen.df)

print(paste0(gcb, "..."), quote=FALSE)
toExport <- c("ref.v", "query.v", "query.v.len", "Cp.v", "ct.v", "chr.v.len",
              "totbin", "gcb", "persist.dir", "out.dir")

#### PARALLEL EXECUTION #########
OLAPCP.MX <- foreach(itr=isplitVector(1:ct.v.len, chunks=nCPU), 
                     .inorder=TRUE, .combine="rbind",
                     .export=toExport, 
                     .noexport=ls()[!ls()%in%toExport]
                     
) %op% {
  
  chunk <- sapply(X=itr, simplify=FALSE, FUN=function(i){
    
    ct <- ct.v[i]
    # Collect bins per chromosome separating per Cp
    binPCP.lst <- vector(mode="list", length=length(Cp.v))
    names(binPCP.lst) <- Cp.v
    
    for(c in 1:chr.v.len){
      
      chr <- chr.v[c]
      chr.num <- gsub(x=chr, pattern="chr", replacement="", fixed=TRUE)
      
      load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
      if(ct%in%Ct.v){
        ct.TF <- PERSIST.MX$hits[[ct]]>0L
      } else if(ct=="allCT"){
        ct.TF <- rep(TRUE, times=length(PERSIST.MX$ntis))
      } 
      cp.v <- unique(PERSIST.MX$ntis[ct.TF])
      for(cp in cp.v){
        x  <- unique(unlist( PERSIST.MX$hits[ct.TF & PERSIST.MX$ntis==cp,c("i", "j")] ))
        cp <- as.character(cp)
        binPCP.lst[[cp]] <- unique( c(binPCP.lst[[cp]], paste(chr.num, x, sep="_")) )
        rm(x)
      }
      rm(PERSIST.MX, cp.v, chr.num, ct.TF); gc()
      
      print(paste0(chr, " data added!"), quote=FALSE)
      
    } # chr.v.len for loop end
    
    binPCP.lst[["HiCAll"]] <- unique(unlist(binPCP.lst))
    bincountPCP <- unlist(lapply(X=binPCP.lst, FUN=length), use.names=TRUE)
    names(bincountPCP) <- paste0(names(bincountPCP), "_denom")
    countOlap.v <- c(HiCAll=length(binPCP.lst[["HiCAll"]]),
                     rep(NA, times=query.v.len))
    names(countOlap.v) <- c("HiCAll", paste0(query.v, "_", ref.v))
    for(k in 1:query.v.len){
      query <- query.v[k]
      ref <- ref.v[k]
      olap <- intersect( binPCP.lst[[ref]], binPCP.lst[[query]] )
      countOlap.v[paste0(query, "_", ref)] <- length(olap)
      rm(olap, query, ref)
    } # olap.v for loop end
    
    rm(binPCP.lst); gc()
    countOlap.v <- c(countOlap.v, `Allbin_denom`=totbin, 
                     bincountPCP[paste0(ref.v, "_denom")])
    if(length(countOlap.v)%%2!=0){ stop("Length not even.") }
    
    print(paste0(ct, " done!"), quote=FALSE)
    return(countOlap.v)
    
  })
  return(do.call("rbind",  chunk))
}
### END OF PARALLEL EXECUTION ###

dimnames(OLAPCP.MX)[[1]] <- ct.v
collen <- ncol(OLAPCP.MX)
denom.mx <- OLAPCP.MX[,(collen/2+1):collen]
denom.mx <- cbind(denom.mx[,"Allbin_denom"], denom.mx)
# Counts
OLAPCP.MX <- cbind(`NonHiC_Allbin`=totbin-OLAPCP.MX[,"HiCAll"], 
                   `HiCAll_Allbin`=OLAPCP.MX[,"HiCAll"],
                   OLAPCP.MX[,2:(collen/2)])
write.csv(x=OLAPCP.MX, file=paste0(out.dir, "/", gcb, "_chrALL_count_ubinsOlapCp.csv"),
          row.names=TRUE, quote=FALSE) 
save(OLAPCP.MX, file=paste0(out.dir, "/", gcb, "_chrALL_count_ubinsOlapCp.RData"))
# Fraction
OLAPCP.MX <- OLAPCP.MX/denom.mx; rm(denom.mx)
write.csv(x=OLAPCP.MX, file=paste0(out.dir, "/", gcb, "_chrALL_fraction_ubinsOlapCp.csv"),
          row.names=TRUE, quote=FALSE) 
save(OLAPCP.MX, file=paste0(out.dir, "/", gcb, "_chrALL_fraction_ubinsOlapCp.RData"))
rm(OLAPCP.MX); gc()  

# rm(list=ls()); gc()




