################################################################################
# In a csv file, list possible ways to group Phylo-HMRF states to be used for 
# plotting
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=T) # Warning for left to right partial matching by $
options(warn=1) # Expands warnings

wk.dir = paste0("/Users/ltamon/SahakyanLab/GenomicContactDynamics/11_Complementarity")
src.file = paste0(wk.dir, "/Phylo-HMRF_state_group_src.csv")
waysToGroup = 3
### OTHER SETTINGS #############################################################
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
df <- read.csv(src.file, header=F, stringsAsFactors=F)
df$V3 <- NULL
add.mx <- matrix(data=NA_character_, nrow=length(df$V1), ncol=waysToGroup) 

grp.num = 1
add.mx[,grp.num] <- df$V2 
for( species in c("hom", "pan", "pyg", "gor") ){
  
  merge.TF <- grepl(x=df$V2, pattern=paste0("NC-", species))
  add.mx[merge.TF,grp.num] <- paste0("NC-", species)
  
}

grp.num = 2
merge.TF <- grepl(x=df$V2, pattern="NC-pan|NC-pyg|NC-gor")
add.mx[,grp.num] <- df$V2 
add.mx[merge.TF,grp.num] <- "NC-1nonhom"

merge.TF <- grepl(x=df$V2, pattern="NC-hom", fixed=T)
add.mx[merge.TF,grp.num] <- "NC-hom"

grp.num = 3
merge.TF <- grepl(x=df$V2, pattern="NC")
add.mx[,grp.num] <- df$V2
add.mx[merge.TF,grp.num] <- "NC"

merge.TF <- grepl(x=df$V2, pattern="C-high|C-mid|C-low")
add.mx[merge.TF,grp.num] <- "C"

#

state.groupings <- cbind(df, add.mx)
colnames(state.groupings) <- c( "State", paste0("Grp", 0:(waysToGroup)) 
                              )
out.file <- gsub(src.file, pattern="_src", replacement="", fixed=T)
write.csv(state.groupings, file=out.file, row.names=F)

# Colours
df <- read.csv(src.file, header=F, stringsAsFactors=F)
df$V1 <- NULL
group.colours <- df[!duplicated(df),]
group.colours <- rbind(group.colours,
                       c("NC-hom", "80B1D2"),
                       c("NC-pan", "E4AB23"),
                       c("NC-pyg", "76AC75"),
                       c("NC-gor", "E289B9"),
                       c("NC-1nonhom", "D7D8D7"),
                       c("NC-hom", "777877"),
                       c("C", "670912")
                       )

colnames(group.colours) <- c("Grp", "Col")

group.colours$Col[group.colours$Grp == "NC"] <- "414241" 
  
out.file <- gsub(src.file, pattern="src", replacement="col", fixed=T)
write.csv(group.colours, file=out.file, row.names=F)

# rm(list=ls()); gc()