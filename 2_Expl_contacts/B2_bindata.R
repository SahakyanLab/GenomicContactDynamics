################################################################################
# Comprehensive bin data information indicating if the bin forms a contact 
# per Cp per tissue. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/19_Circos"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Circos"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
out.dir = paste0(persist.dir, "/binmx")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb.v = c("min2Mb", "min05Mb")
chr.v = paste0("chr", c(1:22, "X"), sep="")
CT.v <- c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB", "AG",
          "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
CP.v = 1:21 
bin.len = 40000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Chromosome length file
chrLen.df <- read.table(file=chrLenfile, as.is=FALSE, header=TRUE,
                        colClasses=c("character", "integer", "integer"))

for(gcb in gcb.v){
  for(chr in chr.v){
    # Load PERSIST.MX
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    ubins <- unique( c(unique(PERSIST.MX$hits$i), unique(PERSIST.MX$hits$j)) )
    ct.v <- sort( colnames(PERSIST.MX$hits)[-(1:2)], decreasing=FALSE )
    cp.v <- sort( unique(PERSIST.MX$ntis), decreasing=FALSE )
    bin.last <- ceiling( chrLen.df[chrLen.df$chromosome==chr,"length.bp"]/bin.len )
    # Initialise matrix
    BIN.MX <- matrix(data=0, nrow=bin.last, ncol=2+length(CT.v)+length(CP.v), 
                     dimnames=list(1:bin.last, c("start", "end", CT.v, CP.v))
    )
    BIN.MX[,"end"] <- (1:bin.last)*bin.len
    BIN.MX[,"start"] <- (BIN.MX[,"end"])-bin.len+1
    BIN.MX[bin.last,"end"] <- chrLen.df[chrLen.df$chromosome==chr,"length.bp"]
    #for(bin in ubins){
    #  # Cp assignment of bins
    #  ij.TF <- PERSIST.MX$hits$i==bin | PERSIST.MX$hits$j==bin
    #  BIN.MX[bin,as.character(1:21)] <- as.integer(CP.v%in%unique(PERSIST.MX$ntis[ij.TF]))
    #  # Ct asignment of bins
    #  ct.TF <- colSums(PERSIST.MX$hits[ij.TF,-(1:2)]) > 0
    #  BIN.MX[bin,CT.v] <- as.integer(ct.TF[CT.v])
    #  rm(ij.TF, ct.TF)
    #}
    for(ct in ct.v){
      bin.incl <- PERSIST.MX$hits[PERSIST.MX$hits[[ct]]>0, c("i", "j")]
      bin.incl <- unique( c(unique(bin.incl$i), unique(bin.incl$j)) )
      BIN.MX[as.numeric(bin.incl),ct] <- 1
    }
    for(cp in cp.v){
      bin.incl <- PERSIST.MX$hits[PERSIST.MX$ntis==cp, c("i", "j")]
      bin.incl <- unique( c(unique(bin.incl$i), unique(bin.incl$j)) )
      BIN.MX[as.numeric(bin.incl),as.character(cp)] <- 1
    }
    rm(PERSIST.MX);gc()
    save(BIN.MX, file=paste0(out.dir, "/", chr, "_", gcb, "_bindata.RData"))
    print(paste0(chr, " done!"), quote=FALSE)
  }
}

# rm(list=ls()) 


