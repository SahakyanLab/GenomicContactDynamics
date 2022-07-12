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
hub.dir = paste0(wk.dir, "/out_hubfile/All_LTr")
out.dir = paste0(wk.dir, "/out_hubsummary/All_LTr")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste("chr", c(1:22, "X"), sep="")
gap.bin.v = 50
hub.id = "All_topCP3_gapBin50"
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
HUB <- list.files(hub.dir)
HUB <- HUB[grepl(x=HUB, pattern=hub.id)]
HUB <- unique(gsub(x=HUB, pattern=".csv|_genes.txt", replacement=""))

HUBSUM <- list()
GENES <- list()
for(chr in chr.v){
  
  chr.TF <- grepl(x=HUB, pattern=paste0("_", chr, "_"), fixed=TRUE)
  
  for(gap.bin in gap.bin.v){
    
    gap.TF <- grepl(x=HUB, pattern=paste0("_gapBin", gap.bin, "_"), fixed=TRUE)
    hub.v <- sort(HUB[chr.TF & gap.TF], decreasing=FALSE)
    
    for(hub in hub.v){
      
      hub.num <- strsplit(x=hub, split="_hub")[[1]][2]
      hub.tbl <- read.csv(file=paste0(hub.dir, "/", hub, ".csv"), header=TRUE, 
                          stringsAsFactors=FALSE)
      Nij.center <- max(hub.tbl$Nij)
      hub.center <- hub.tbl[hub.tbl$Nij==Nij.center, "bin"]
      genes <- hub.tbl$gene
      
      if( all(is.na(genes)) ){
        
        genes <- NA
        Ngenes <- 0
        print(paste0("No genes for ", hub), quote=FALSE)
        
      } else {
        
        genes <- unlist(strsplit(x=genes, split=";")) 
        genes <- unique(genes[!is.na(genes)])
        Ngenes <- length(genes)
        
      }
      # Save genelist for functional annotation
      writeLines(genes, con=paste0(out.dir, "/", hub, "_genes.txt"))
      GENES[[hub]] <- genes
      
      HUBSUM[[hub]] <- c(hub=hub, chr=chr, gap=gap.bin, hubN=hub.num, center=hub.center, 
                         NijCenter=Nij.center, Ngenes=Ngenes, genes=paste(genes, collapse=";"))
      
      print(paste0(hub, " done!"), quote=FALSE)
      
    }
  }
}

NGENES <- unlist(lapply(X=GENES, FUN=length))
GENES.MX <- matrix(data=NA, nrow=max(NGENES), ncol=length(GENES),
                   dimnames=list(NULL, names(GENES)))
for(hub in names(GENES)){
  GENES.MX[1:NGENES[hub],hub] <- GENES[[hub]]
}
write.table(GENES.MX, file=paste0(out.dir, "/", gcb, "_", hub.id, "_hubgenes.txt"),
            sep="\t", quote=FALSE)

HUBSUM <- do.call("rbind.data.frame", c(HUBSUM, stringsAsFactors=FALSE))
colnames(HUBSUM) <- c("hub", "chr", "gap", "hubN", "center", "NijCenter", 
                      "Ngenes", "genes")
write.csv(HUBSUM, file=paste0(out.dir, "/", gcb, "_", hub.id, "_hubsum.csv"),
          row.names=FALSE, quote=FALSE)

# rm(list=ls())



