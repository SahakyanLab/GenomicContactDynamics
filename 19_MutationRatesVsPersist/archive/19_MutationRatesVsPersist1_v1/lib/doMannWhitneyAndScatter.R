############################################################################### 
# Specific function to do Mann-Whitney test using x and df object.
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
#library(ggplot2)
#source(paste0(lib, "/GG_bgr.R"))
#source(paste0(wk.dir, "/lib/identifyAltHyp.R"))
### FUNCTION ###################################################################
doMannWhitneyAndScatter <- function(x=x, df=df, calc=calc, plot.id=id
                                    ){
  
  # Wilcoxon-Mann-Whitney test (non-parametric, unpaired/independent groups) 
  # to compare orig and shuff values per Cp
  # The test assumes that the shape of the two distributions are similar so
  # check the boxplots
  # If both x and y are given and paired is FALSE, a Wilcoxon rank sum test
  # (equivalent to the Mann-Whitney test) is carried out.
  df$pval <- df$alt <- NA
  ind.v <- sort(unique(df$ind), decreasing=FALSE)
  
  if( typeof(x$ind)!=typeof(df$ind) ){
    stop("x$ind and df$ind of different type.")
  }
  
  for(ind in ind.v){
    
    id.mean <- paste0(ind, "_MEAN")
    id.med <- paste0(ind, "_MEDIAN")
    
    alt.v <- identifyAltHyp(x=c(df[id.mean,"values"], df[id.med,"values"]),
                            y=c(df["1_MEAN","values"], df["1_MEDIAN","values"])
                            )
    
    # The one-sided alternative "greater" is that x is shifted to the right of y
    df[id.mean,"pval"] <- wilcox.test(x=x[x$ind==ind,calc], y=x[x$ind==1,calc],
                                      alternative=alt.v["MEAN"], paired=FALSE)$p.value 
    df[id.med,"pval"] <- wilcox.test(x=x[x$ind==ind,calc], y=x[x$ind==1,calc],
                                     alternative=alt.v["MEDIAN"], paired=FALSE)$p.value 
    
    df[id.mean,"alt"] <- alt.v["MEAN"] 
    df[id.med,"alt"] <- alt.v["MEDIAN"] 
    
    rm(alt.v, id.mean, id.med)
    
  }
  
  # Plot
  df$pval <- c(`0`="ns", `1`="sig")[ as.character(as.numeric(df$pval<0.05)) ]
  df$pval <- factor(df$pval, levels=c("ns", "sig"))
  df$ind <- factor(df$ind, levels=c(as.character(sort(unique(df$ind)))))
  
  df$alt <- factor(as.character(df$alt), levels=c("greater", "less", "two.sided"))
  df$id <- factor(as.character(df$id), levels=c("MEAN", "MEDIAN"))
  
  p <- ggplot(data=df, aes(x=ind, y=values, group=id)) +
    geom_point(aes(colour=id, alpha=pval, shape=alt), size=3) +
    scale_colour_manual(values=c(MEAN="#F8766D", MEDIAN="#00BFC4")) + 
    scale_alpha_manual(values=c(ns=0.3, sig=1)) + 
    scale_shape_manual(values=c(greater=16, less=17, `two.sided`=15)) + 
    labs(x="Cp", y=calc, title=plot.id) + 
    bgr2 + 
    theme(plot.title=element_text(size=2))
  
  return(list(df=df, p=p))
  
}
