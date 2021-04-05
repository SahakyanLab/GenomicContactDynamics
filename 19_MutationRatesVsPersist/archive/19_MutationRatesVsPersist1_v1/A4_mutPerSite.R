################################################################################
# Plot:
# a. Total mutations per ID_SAMPLE in the dataset. Also calculates (i. )the mean 
# and median of the total mutations per sample and (ii.) the number of unique 
# ID_SAMPLE.
# b. Total mutations per nucleotide site per chromosome. Make plot for all and 
# each mutation type (e.g. A>C). 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster"  # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

# Expands warnings
options(warn=1)

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/t1-data/user/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/19_Mutation_rates"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
#src.dir = paste0(wk.dir, "/out_filter")
src.dir = paste0(data.dir, "/cancer_somatic_mutation/out_dataPerSignature")
out.dir = paste0(wk.dir, "/out_mutPerSite")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
data.id = "donor_centric_PCAWG_as_is" # "CosmicNCV", "donor_centric_PCAWG_sigEperc50_MMR12"
src.id = "Hg19" # "hg38ToHg19" | "Hg19" 
chr.v = NULL
bin.len = 40000
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(ggsci)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
chrlen.df <- read.delim(file=chrLenfile, header=TRUE)

load(file=paste0(src.dir, "/", data.id, "_", src.id, "_final.RData"))
ncv.df <- ncv.df[ncv.df$Coding==0,] # Non-coding mutations only
ncv.df <- ncv.df[,c("ID_SAMPLE", "chr", "start", "MUT")]

# Plot total mutation (plus mean and median) per ID_SAMPLE
ID_SAMPLE.len <- length(unique(ncv.df$ID_SAMPLE))

x <- table(ncv.df$ID_SAMPLE)
med <- median(x)
men <- mean(x)

x <- data.frame(ID_SAMPLE=as.character(ncv.df$ID_SAMPLE), stringsAsFactors=FALSE)
p <- ggplot(data=x, aes(x=ID_SAMPLE)) + 
  geom_bar() +
  geom_hline(yintercept=med, linetype="dotted", color="red", size=2) +
  geom_hline(yintercept=men, linetype="dotted", color="blue", size=2) +
  labs(x="ID_SAMPLE", y="# Mut", 
       title=paste0(data.id, "_", src.id, "_", ID_SAMPLE.len, "uniqID_SAMPLE_med(red)=", med,
                    "_mean(blue)=", men, "\nNfinalnoncodingmut=", nrow(ncv.df))
       ) 
ggsave(filename=paste0(out.dir, "/", data.id, "_", src.id, "_final_mutPerSamp.pdf"),
       height=10, width=10, plot=p)
rm(ID_SAMPLE.len, med, men, p, x); gc()

#-------------------
if( is.null(chr.v) ){ chr.v <- unique(ncv.df$chr) }
for(chr in chr.v){
  
  chr.len <- chrlen.df$length.bp[chrlen.df$chromosome==chr]
  
  incl.TF <- ncv.df$chr==chr
  
  p <- ggplot(data=ncv.df[incl.TF, c("start", "MUT")], aes(x=start)) +
    geom_bar(aes(fill=MUT)) +
    scale_x_continuous(limits=c(1,chr.len)) +
    scale_fill_npg() + 
    labs(x="Position", y="# Mut", fill="Base", title=paste0(data.id, "_", src.id, "_", chr)) 
  
  ggsave(filename=paste0(out.dir, "/", chr, "_", data.id, "_", src.id, "_final_mutPerSite.pdf"),
         height=10, width=20, plot=p)
  
  # Panel, per mutation type
  p <- ggplot(data=ncv.df[incl.TF, c("start", "MUT")], aes(x=start)) +
    geom_bar(aes(fill=MUT)) +
    scale_x_continuous(limits=c(1,chr.len)) +
    scale_fill_npg() + 
    labs(x="Position", y="# Mut", fill="Base", title=paste0(data.id, "_", src.id, "_", chr)) +
    facet_grid(MUT~.)
  
  ggsave(filename=paste0(out.dir, "/", chr, "_", data.id, "_", src.id, "_final_panel_mutPerSite.pdf"),
         height=10, width=20, plot=p)
  
  ncv.df <- ncv.df[!incl.TF,]
  rm(p, chr.len); gc()
  print(paste0(chr, " done!"), quote=FALSE)
  
}

# rm(list=ls()); gc()

