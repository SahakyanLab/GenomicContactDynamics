################################################################################
### FUNCTION ###################################################################
applyFilter3 <- function(ncv.df=ncv.df){
  
  warning("Make sure ncv.df was filtered for only Whole_Genome_Reseq or Whole_Exome 
          mutations and was applied with applyFilter1() and applyFilter2.")
  
  x <- list()
  
  # zygosity?
  # Removing SNP filter only adds ~60K entries
  #ncv.df$SNP=="n"
  
  base.v <- c("A", "T", "C", "G")
  base.v <- c(base.v, tolower(base.v))
  
  x$incl.TF <- ncv.df$`Mutation somatic status`=="Confirmed somatic variant" &
               ncv.df$WT_SEQ%in%base.v & ncv.df$MUT_SEQ%in%base.v &
               (ncv.df$end-ncv.df$start)==0 & 
               ncv.df$chr%in%paste0("chr", c(1:22)) &
               ncv.df$Whole_Genome_Reseq=="y" & ncv.df$Whole_Exome=="n"
  
  return(x)

}

# rm(list=ls()); gc()

