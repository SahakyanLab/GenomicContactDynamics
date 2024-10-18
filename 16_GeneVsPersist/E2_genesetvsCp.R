################################################################################
# Check overlap of housekeeping gene list with gene working tables
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
lib = file.path(home.dir, "git/GenomicContactDynamics/lib/")
data.dir = file.path(home.dir, "data/gcd/Database")
annofile.dir = file.path(data.dir, "ucsc_tables/hsa_geneAnno")
wk.dir = file.path(home.dir, "git/GenomicContactDynamics/16_GeneVsPersist")
annofilePath = file.path(annofile.dir, "hg19anno_ALL")
# Use longest transcript to remove bias from uneven count of transcripts per gene
#annofilePath = paste0(annofile.dir, "/hg19annoLTr_ALL")
# Unique HUGO genes overlapping with Cp regions (any type of overlap)
CpGenesPath = file.path(home.dir, "data/gcd/16_GeneVsPersist/out_anno_union/min2Mb_ALL_name2")
hkgenes.dir = file.path(data.dir, "housekeeping_genes/human")
out.dir = file.path(wk.dir, "out_genesetvsCp")
### OTHER SETTINGS #############################################################
# Annotation file prefix
anno.nme = "hg19anno" #"hg19annoLTr"
gcb="min2Mb"
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(ggplot2)
source(file.path(lib, "do_permute_test.R"))
source(file.path(lib, "GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
anno.df <- read.delim(file = annofilePath, header = TRUE, stringsAsFactors = FALSE)
load(file.path(hkgenes.dir, "Housekeeping_GenesHuman.RData"))
hkgenes_df <- Housekeeping_Genes
hkgenes_top_df <- read.delim(file.path(hkgenes.dir, "MostStable.csv"), sep = ";")
hkgenes_df$top <- ifelse(hkgenes_df$Gene.name %in% hkgenes_top_df$Gene.name, "yes", "no")

# Housekeeping genes
genes_of_interest <- unique(as.character(hkgenes_df$Gene.name))
#genes_of_interest <- unique(as.character(hkgenes_df$Gene.name[hkgenes_df$top == "yes"]))

# Retrieve Cp genes
temp <- readLines(con=CpGenesPath)
CpGenesInd <- grep(x=temp, pattern=">all_genes_cp_", fixed=TRUE)
CpGenes <- temp[CpGenesInd + 1L]
names(CpGenes) <- gsub(x=temp[CpGenesInd], pattern=">all_genes_cp_|_end", replacement="")
rm(temp, CpGenesInd); gc()
CpGenes.len <- length(CpGenes)
# List of genes overlapping with all long-range contact regions and regions
# of Cp 1 to 21. 
CpGenes <- strsplit(x=CpGenes, split=";", fixed=TRUE)

obs_df <- sapply(X=names(CpGenes), simplify=FALSE, FUN=function(cp){
  if( any(duplicated(CpGenes[[cp]])) ){ stop(paste0("Duplicated genes in ", cp)) }
  n0 <- length(CpGenes[[cp]])
  frgoi <- round(sum(CpGenes[[cp]] %in% genes_of_interest)/n0, digits=3)
  frrest <- round(sum(!CpGenes[[cp]] %in% genes_of_interest)/n0, digits=3)
  cp <- as.character(cp)
  if( (frgoi+frrest)==1 ){
    print(paste0(cp, " done!"), quote=FALSE)
    return( cbind.data.frame(Cp=c(cp,cp), type=c("rest", "goi"), fr=c(frrest, frgoi)) )
  } else {
    stop("Sum is not 1.")
  }
})
obs_df <- do.call("rbind.data.frame", obs_df)

# Convert Cp column to factor, arranging the levels for the plot
obs_df$Cp <- as.character(obs_df$Cp)
obs_df$Cp <- factor(x=obs_df$Cp, levels=unique(obs_df$Cp))

p <- ggplot(data=obs_df, aes(x=Cp, y=fr, fill=type)) +
  geom_col() +
  labs(title=paste0(gcb, "_", anno.nme), x=expression(bold("c"["p"])),
       y="Fraction", fill=NULL) + 
  bgr2 +
  theme(axis.text.x = element_text(size=10, angle=360, colour="black"))
ggsave(filename=paste0(out.dir, "/", gcb, "_", anno.nme, "_genesetvsCp.pdf"),
       plot = p, width=10, height=10)

# Permutation test to determine significance of fraction of genes of interest per Cp
p_values <- numeric(length(CpGenes))

p_values <- vapply(names(CpGenes), FUN.VALUE = numeric(1), FUN = function(x) {
  do_permute_test(unique(CpGenes$HiC_all), 
                  obs_set = unique(unlist(CpGenes[as.character(x)])), 
                  olap_set = unique(genes_of_interest),
                  seed_val = 290, n_iter = 1000)
})
p_adj <- p.adjust(p_values, method = "BH")

# Scatter plot

plot_df <- obs_df %>% 
  filter(type == "goi")

if (identical(names(p_adj), as.character(plot_df$Cp))) {
  plot_df <- obs_df %>% 
    filter(type == "goi") %>% 
    mutate(p_adj = p_adj, perc = fr * 100)
}

p <- ggplot(plot_df, aes(x = Cp, y = perc)) +
  geom_point(aes(colour = p_adj < 0.05)) +
  scale_colour_manual(values = c("black", "darkred")) +
  labs(y = "Percentage overlap", title = "Overlap of housekeeping genes") + 
  bgr1 +
  theme(axis.text.x = element_text(size = 5))

ggsave(filename=paste0(out.dir, "/", gcb, "_", anno.nme, "_genesetvsCp_fraction.pdf"),
       plot = p, width=5, height=5)

# ggsave(filename=paste0(out.dir, "/", gcb, "_", anno.nme, "_genesetvsCp_fraction_tophousekeeping.pdf"),
#        plot = p, width=5, height=5)

# rm(list=ls()); gc()

