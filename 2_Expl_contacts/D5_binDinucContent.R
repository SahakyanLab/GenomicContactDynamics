################################################################################
# Plot GC/AT content vs. Cp. Deliberately excluded last bin to keep bin length to
# 40 kb to keep working with base counts instead of fraction until needed.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/SahakyanLab/GenomicContactDynamics/2_Expl_contacts"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
persist.dir = paste0(data.dir, "/HiC_features_GSE87112_RAWpc")
binkmer1.dir = paste0(persist.dir, "/out_binBaseContent/maskfile0")
out.dir  = paste0(wk.dir, "/out_binDinucContent")
chrlen.file = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = paste0("chr", c(1:22, "X")) 
bin.len = 40000
plotOnly = T
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if(plotOnly==F){
  
  chrlen.df <- read.delim(chrlen.file, stringsAsFactors=F, header=T)
  
  for(chr in chr.v){
    
    chr.len <- chrlen.df$length.bp[chrlen.df$chromosome==chr]
    
    load(paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    
    binPerCp <- by(data=PERSIST.MX$hits[,c("i", "j")],
                   INDICES=PERSIST.MX$ntis,
                   FUN=function(bin.df) unique(unlist(bin.df)))
    
    load(paste0(binkmer1.dir, "/", chr, "_BinKmer1.RData"))
    
    # Exclude last bin with length < HiC res, all with missing base anyway except
    # for chr17
    last.bin <- ceiling(chr.len/bin.len)
    if( identical(as.numeric(1:last.bin), as.numeric(BINKMER.MX[,"bins"])) ){
      BINKMER.MX <- BINKMER.MX[-(last.bin),] 
    } else {
      stop(chr, ": Missing bins in BINKMER.MX")
    }
    
    tmp <- unique(na.omit(rowSums(BINKMER.MX[,c("A", "C", "G", "T")], na.rm=F)))
    
    if( !identical(tmp, bin.len) ){
      stop(chr, ": Faulty ATGC content.")
    } else {
      
      GCs <- rowSums(BINKMER.MX[,c("C","G")], na.rm=F)
      rm(BINKMER.MX)
      
      GCPerCp <- sapply(X=names(binPerCp), simplify=F, FUN=function(Cp.nme){
        
        bins <- binPerCp[[Cp.nme]]
        GC <- GCs[bins]
        return(GC)
        
      })
      
      rm(binPerCp, GCs)
      
      save(GCPerCp, file=paste0(out.dir, "/", chr, "_GC_tmp_", gcb, ".RData"))
      
      rm(GCPerCp)
      gc()
      
    }
    
    print(paste0(chr, " done!"), quote=F)
    
  }
  
} else {
  
  GC <- setNames(object=vector("list", 21), nm=1:21)
  
  for(chr in chr.v){
    
    tmp <- setNames(object=vector("list", 21), nm=1:21)
    
    load(file=paste0(out.dir, "/", chr, "_GC_tmp_", gcb, ".RData"))
    
    tmp[names(GCPerCp)] <- GCPerCp 
    
    GC <- Map(c, GC[names(GC)], tmp[names(GC)])
    
    rm(tmp, GCPerCp)
    
  }
  
  print(lengths(GC), quote=F)
  
  GC <- stack(GC)
  AT <- GC
  AT$values <- bin.len - GC$values
  
  out.name <- paste0(gcb, "_chrALL_dinucContentVsCp")
  png(file=paste0(out.dir, "/", out.name, ".png"), height=300*5, width=300*10, res=300)
  par(mfrow=c(1,2))
  
  for( dinuc in c("GC", "AT") ){
    
    eval(parse(text=paste0("df <- ", dinuc)))
    
    df$values <- df$values/bin.len
    boxplot(formula=values ~ ind, data=df, outline=T, xlab="Cp", ylab="Fr", col="#FDC776", cex=0.2,
            main=paste0(out.name, "_", dinuc, "_lastbinexcluded"), cex.main=0.5, ylim=c(0.2,0.8))
    
  }
  
  dev.off()
  
}

# rm(list=ls()); gc()