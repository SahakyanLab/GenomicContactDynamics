############################################################################### 
# Specific function to do Mann-Whitney test using x and df object.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#source(paste0(wk.dir, "/lib/identifyAltHyp.R"))
### FUNCTION ###################################################################
doMannWhitney <- function(x=x, df=df, calc=calc){
  
  # Wilcoxon-Mann-Whitney test (non-parametric, unpaired/independent groups) 
  # to compare orig and shuff values per Cp
  # The test assumes that the shape of the two distributions are similar so
  # check the boxplots
  # If both x and y are given and paired is FALSE, a Wilcoxon rank sum test
  # (equivalent to the Mann-Whitney test) is carried out.
  df$pval <- df$alt <- NA
  
  if( !is.character(x$ind) | !is.character(df$ind) ){
    stop("x$ind and/or df$ind not character as expected.")
  }
  
  ind.v <- as.character(unique(df$ind))
  rownames(df) <- as.character(df$ind)
  
  for(ind in ind.v){
    
    row.id <- as.character(ind)
    alt <- unname(identifyAltHyp(x=df[row.id,"values"], y=df["1","values"]))

    # The one-sided alternative "greater" is that x is shifted to the right of y
    df[row.id,"pval"] <- wilcox.test(x=x[x$ind==ind,calc], y=x[x$ind=="1",calc],
                                     alternative=alt, paired=FALSE)$p.value 
    df[row.id,"alt"] <- alt 
   
  }
 
  return(df=df)
  
}
