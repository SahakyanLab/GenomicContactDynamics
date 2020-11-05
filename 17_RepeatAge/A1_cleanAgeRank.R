################################################################################
# Compare Giordano_405HumanTEs and Published_136 data age ranks
# Get repeat subfamilies from age ranks that overlap with Repeat Masker
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"
if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/17_RepeatAge"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(wk.dir, "/out_cleanAgeRank")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# All list have unique repNames
Giordano <- fread(file=paste0(data.dir, "/Giordano2007_405HumanTEs_renamed.txt"),
                     header=FALSE, data.table=FALSE, stringsAsFactors=FALSE)[[1]]
Published.df <- fread(file=paste0(data.dir, "/publishedAges_136sorted.csv"), sep=",",
                         header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)
Kapusta <- fread(file=paste0(data.dir, "/Kapusta2017_537TEsandFamilyorigin_unique529.csv"),
                 header=TRUE, data.table=FALSE, sep=",", stringsAsFactors=FALSE)
repmasker <- unique( fread(file=paste0(data.dir, "/hg19Rep_classfamsubfam_unique.txt"),
                           header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)[["repName"]]
                     )

# Remove repNames not in repmasker
# Order kept
#----------------------------------------
Giordano.get <- intersect(Giordano, repmasker)
write(x=Giordano.get, file=paste0(data.dir, "/Giordano2007_", length(Giordano.get), 
                              "HumanTEs_renamed.txt"))
Giordano.rm <- setdiff(Giordano, Giordano.get)
write(x=Giordano.rm, file=paste0(data.dir, "/Giordano2007_", length(Giordano.rm), 
                                  "HumanTEs_renamed_REMOVED.txt"))
#----------------------------------------
Published.df <- Published.df[order(Published.df$ageMyrConsensus, decreasing=TRUE),]
Published.get.df <- Published.df[Published.df$repNameRenamed %in% repmasker, ]
write(x=Published.get.df$repNameRenamed, 
      file=paste0(data.dir, "/Published_", nrow(Published.get.df), 
                  "HumanTEs_renamed.txt"))
write.csv(x=Published.get.df, file=paste0(data.dir, "/Published_", nrow(Published.get.df), 
                                    "HumanTEs_renamed.csv"),
          row.names=FALSE, quote=FALSE)
Published.rm <- setdiff(Published.df$repNameRenamed, Published.get.df$repNameRenamed)
write(x=Published.rm, 
      file=paste0(data.dir, "/Published_", length(Published.rm),
                  "HumanTEs_renamed_REMOVED.txt"))

orderPbasedonG <- match(Published.get.df$repNameRenamed, Giordano.get)
mx <- cbind(repName = Published.get.df$repNameRenamed, orderOnGiordano=orderPbasedonG )
# OrderPbasedonG <- orderPbasedonG[!is.na(orderPbasedonG)]
write.table(mx, file=paste0(data.dir, "/orderPublished_basedonGiordano.txt"),
            col.names=TRUE, row.names=FALSE, quote=FALSE)
#----------------------------------------
# Overlap Repeat masker and Kapusta2017 data
Kapusta.get <- Kapusta[!is.na( match(Kapusta$Repeat_Name, repmasker) ),]
write.table(Kapusta.get, file=paste0(data.dir, "/Kapusta2017_", nrow(Kapusta.get), 
                                 "TEsandFamilyorigin_unique.txt"),
            col.names=TRUE, row.names=FALSE, quote=FALSE)
Kapusta.rm <- setdiff(Kapusta$Repeat_Name, Kapusta.get$Repeat_Name)
write(Kapusta.rm, file=paste0(data.dir, "/Kapusta2017_", length(Kapusta.rm), 
                              "TEsandFamilyorigin_unique_REMOVED.txt"))
#----------------------------------------
# Add published ages and Giordano360retro ranking to GiorPubl ranking
Published.get.df <- fread(file=paste0(data.dir, "/Published_105HumanTEs_renamed.csv"),
                          header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)
GiorPubl372 <- fread(file=paste0(data.dir, "/GiordanoANDPublished_372HumanTEs_renamed.txt"),
                     header=FALSE, data.table=FALSE, stringsAsFactors=FALSE)[[1]]
ind1 <- match(x=GiorPubl372, table=Published.get.df$repNameRenamed)

# Giordano et al 2007 ranking for 360 retrotransposons
Giordano360RetroRank <- fread(file=paste0(data.dir, "/Giordano2007_360HumanTEs.csv"),
                           header=TRUE, data.table=FALSE, stringsAsFactors=FALSE)
ind2 <- match(x=GiorPubl372, table=Giordano360RetroRank$Name)
GiorPubl372withAges <- cbind(repName=GiorPubl372, Published.get.df[ind1,-1],
                             Giordano360RetroRank=Giordano360RetroRank[ind2, "Position"])

write.csv(GiorPubl372withAges, file=paste0(data.dir, "/GiordanoANDPublished_372HumanTEs_renamed_withAges.csv"),
          row.names=FALSE, quote=FALSE)
      
# rm(list=ls()); gc()
