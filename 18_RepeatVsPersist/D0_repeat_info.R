################################################################################
# Dataframe containing repeat subfamilies (repNames) and their repeat family and 
# class. Also contains grouping of repeat subfamilies based on GROUP.CLASS argument.
# There are 1395 unique subfamilies but output contains 1397 rows because there 
# are two repeat subfamilies classified into two families and classes.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon" #"/stopgap/sahakyanlab/" #"/t1-data/user"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/18_RepeatVsPersist")
out.dir = paste0(wk.dir, "/z_ignore_git/out_combine") #paste0(wk.dir, "/out_metric_extremes")
repeatmx.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/17_RepeatAge/z_ignore_git/out_addToSummary")
repeat.file = paste0(repeatmx.dir, "/hg19repeats_repName.RData")
### OTHER SETTINGS #############################################################
GROUP.CLASS <- list(
  not.transposon = c("Low_complexity", "RNA", "rRNA", "Satellite", "scRNA", 
                     "Simple_repeat", "snRNA", "srpRNA", "tRNA"),
  DNA.transposon = c("DNA", "DNA?"), 
  retro.transposon = c("LINE", "SINE", "LTR", "RC", "LINE?", "SINE?", "LTR?"),
  not.classified = c("Other", "Unknown", "Unknown?")
)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
load(file=repeat.file)
rownames(REPEAT.MX) <- NULL
df <- REPEAT.MX[,c("repName", "repClass", "repFamily")]

# Add grouping of repeat class

grp.class.df <- stack(GROUP.CLASS)
grp.class.df$ind <- as.character(grp.class.df$ind)
setnames(grp.class.df , old=c("values", "ind"), new=c("repClass", "plot.group"))

df <- merge(x=df, y=grp.class.df, by="repClass", all.x=T)
rm(grp.class.df)

#

df <- df[order(df$repName, decreasing=F),]
rownames(df) <- NULL

save(df, file=paste0(out.dir, "/repeat_group.RData"))

# rm(list=ls()); gc()