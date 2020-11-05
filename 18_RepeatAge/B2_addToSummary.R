################################################################################
# Add additional ranking, origin, age, reference data to REPEAT.MX
# merge() heavily used here, note that merge does not keep original order of tables
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/3_RepeatAge"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
agerank.dir = paste0(wk.dir, "/out_cleanAgeRank")
repeatmx.dir = paste0(wk.dir, "/out_hg19Repeats_summary")
out.dir = paste0(wk.dir, "/out_addToSummary")
### OTHER SETTINGS #############################################################
# Age rank
Giordano364file = paste0(agerank.dir, "/Giordano2007_364HumanTEs_renamed.txt")
GiorPubl372file = paste0(agerank.dir, "/GiordanoANDPublished_372HumanTEs_renamed.txt")

# published ages file
publAgesfile = paste0(agerank.dir, "/Published_105HumanTEs_renamed.csv")

# Kapusta2017 repeat origin
Kapusta472file = paste0(agerank.dir, "/Kapusta2017_472TEsandFamilyorigin_unique.txt")

# Clusters from RepeatFamiliesVsPersist
# acro not part of clusters 
#ClusterI   <- c("CR1", "hAT", "Deu", "tRNA", "srpRNA", "scRNA")
#ClusterII  <- c("centr", "L1?", "Unknown?")
#ClusterIII <- c("Helitron", "hAT?")
#ClusterIV  <- c("L1", "TcMar-Mariner", "Dong-R4", "TcMar", "TcMar-Tc2", "RTE-BovB",
#                "DNA", "SINE?", "Helitron?", "ERV")
#ClusterV   <- c("L2", "ERVL", "ERVL-MaLR", "MIR", "hAT-Charlie", "ERV1", "RTE",
#                "Low_complexity", "LTR", "TcMar-Tigger", "Gypsy?", "Simple_repeat",
#                "hAT-Tip100", "Gypsy", "Unknown", "DNA?", "TcMar?", "SINE")
#ClusterVI  <- c("Other", "ERVK", "telo", "snRNA", "Satellite", "LTR?", 
#                "Merlin", "Penelope?") #seems to be old?
#ClusterVII <- c("Alu", "ERVL?", "MuDR", "hAT-Blackjack", "rRNA", "PiggyBac",
#                "RNA", "PiggyBac?")
ClusterI = c("L2","L1","ERVL","TcMar-Mariner","ERVL-MaLR","MIR","hAT-Charlie","ERV1",
             "CR1","RTE","Low_complexity","LTR","TcMar-Tigger","Gypsy?","Dong-R4",
             "hAT-Tip100","hAT","Gypsy","Deu","centr","Unknown","ERVL?","MuDR","tRNA",
             "DNA?","hAT-Blackjack","TcMar-Tc2","PiggyBac","DNA","TcMar?","Helitron",
             "hAT?","SINE","SINE?","ERV","LTR?","Penelope?","acro")
ClusterII = c("Alu","Other","Simple_repeat","ERVK","telo","snRNA","Satellite",
              "srpRNA","TcMar","rRNA","RTE-BovB","L1?","scRNA","RNA","PiggyBac?",
              "Helitron?","Merlin","Unknown?")
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(data.table)
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# REPEAT.MX (base)
load(file=paste0(repeatmx.dir, "/hg19repeats_repName_base.RData"))

# add general classification of repeat classes (Kojima 2018 RepBase paper)
REPEAT.MX <- within(data=REPEAT.MX, {
  repType <- NA
  repType[repClass %in% c("SINE", "LINE", "LINE?", "SINE?")] <- "non-LTR"
  repType[repClass=="Other"] <- "non-LTR"
  repType[repClass=="DNA?" & repName %in% c("Eulor6A", "Eulor5A", "Eulor6D", 
                                            "Eulor6C", "Eulor6B", "Eulor5B", 
                                            "Eulor6E")] <- "DNAtransposon"
  repType[repClass %in% c("DNA", "RC")] <- "DNAtransposon"
  repType[repClass=="DNA?" & !(repName %in% c("Eulor6A", "Eulor5A", "Eulor6D", 
                                            "Eulor6C", "Eulor6B", "Eulor5B", 
                                            "Eulor6E"))] <- "Ancient_unclassified"
  repType[repClass %in% c("LTR?", "Unknown", "Unknown?")] <- "Ancient_unclassified"
  repType[repClass=="LTR"] <- "LTR"
  repType[repClass=="Low_complexity"] <- "Low_complexity"
  repType[repClass %in% c("rRNA", "scRNA", "snRNA", "srpRNA", "tRNA", "RNA")] <- "RNA"
  repType[repClass=="Simple_repeat"] <- "Simple_repeat"
  repType[repClass=="Satellite"] <- "Satellite"
})

#----------------------------------
# add clusters from RepeatVsPersist

Clusters.lst <- list(ClusterI=ClusterI,
                     ClusterII=ClusterII
                     #,
                     #ClusterIII=ClusterIII,
                     #ClusterIV=ClusterIV,
                     #ClusterV=ClusterV,
                     #ClusterVI=ClusterVI,
                     #ClusterVII=ClusterVII
                     )

# assign repeat entries to clusters
REPEAT.MX <- within(data=REPEAT.MX, {
  cluster <- NA
  cluster[repFamily %in% ClusterI]   <- "ClusterI"
  cluster[repFamily %in% ClusterII]  <- "ClusterII"
  #cluster[repFamily %in% ClusterIII] <- "ClusterIII"
  #cluster[repFamily %in% ClusterIV]  <- "ClusterIV"
  #cluster[repFamily %in% ClusterV]   <- "ClusterV"
  #cluster[repFamily %in% ClusterVI]  <- "ClusterVI"
  #cluster[repFamily %in% ClusterVII] <- "ClusterVII"
  # cluster=NA:acro repFamily not part of clusters from heatmap (index 318)
  #cluster[repFamily=="acro"] <- "acro"
})

#----------------------------------
# add Giordano et al. (2007) rank order 
# age ranks already ordered: oldest -> youngest

# no duplicates but just in case
Giordano.rank <- unique( fread(file=Giordano364file, header=FALSE, data.table=FALSE, 
                         stringsAsFactors=FALSE)[[1]] )
Giordano.rank <- cbind.data.frame( Giordano364rank=as.integer(1:length(Giordano.rank)), 
                                   repName=Giordano.rank, stringsAsFactors=FALSE)
if( length(setdiff(Giordano.rank$repName, REPEAT.MX$repName))!=0 ){
  stop("Checkpoint 1.")
}

REPEAT.MX <- merge(x=REPEAT.MX, y=Giordano.rank, by="repName", all.x=TRUE) 
REPEAT.MX[is.na(REPEAT.MX[,"Giordano364rank"]),"Giordano364rank"] <- 0L
if( sum(REPEAT.MX[,"Giordano364rank"]!=0) != nrow(Giordano.rank) ){
  stop("Checkpoint 2.")
}

GiorPubl.rank <- unique( fread(file=GiorPubl372file, header=FALSE, data.table=FALSE,
                         stringsAsFactors=FALSE)[[1]] )
GiorPubl.rank <- cbind.data.frame( GiorPubl372rank=as.integer(1:length(GiorPubl.rank)), 
                                   repName=GiorPubl.rank, stringsAsFactors=FALSE)
if( length(setdiff(GiorPubl.rank$repName, REPEAT.MX$repName))!=0 ){
  stop("Checkpoint 3.")
}
REPEAT.MX <- merge(x=REPEAT.MX, y=GiorPubl.rank, by="repName", all.x=TRUE) 
REPEAT.MX[is.na(REPEAT.MX[,"GiorPubl372rank"]),"GiorPubl372rank"] <- 0L
if( sum(REPEAT.MX[,"GiorPubl372rank"]!=0) != nrow(GiorPubl.rank) ){
  stop("Checkpoint 4.")
}

#----------------------------------
# animal origin (eg. eutheria, mammalia etc.)
# 57 repNames are not in Repeat Masker table (see Kapusta 2017 excel file below)
repName_origin <- fread(file=Kapusta472file, header=TRUE, data.table=FALSE,
                        stringsAsFactors=FALSE)
# remove duplicates if present
repName_origin <- repName_origin[!duplicated(x=repName_origin),]

setnames(x=repName_origin, old="Repeat_Age", new="origin")

REPEAT.MX <- merge(x=REPEAT.MX, y=repName_origin[,c("Repeat_Name", "origin")],
                   by.x="repName", by.y="Repeat_Name", all.x=TRUE) 

#----------------------------------
# add published ages gathered
# already ordered
publAges.df <- fread(file=publAgesfile, header=TRUE, data.table=FALSE, 
                     stringsAsFactors=FALSE)
publAges.df <- cbind(Publrank=1:nrow(publAges.df), publAges.df)
if( length(setdiff(publAges.df$repNameRenamed, REPEAT.MX$repName))!=0 ){
  stop("Checkpoint 5.")
}
REPEAT.MX <- merge(x=REPEAT.MX, y=publAges.df, all.x=TRUE, 
                   by.x="repName", by.y="repNameRenamed")
REPEAT.MX[is.na(REPEAT.MX[,c("Publrank")]), "Publrank"] <- 0L
if( sum(REPEAT.MX[,"Publrank"]!=0) != nrow(publAges.df) ){
  stop("Checkpoint 6.")
}
save(REPEAT.MX, file=paste0(out.dir, "/hg19repeats_repName.RData"))
write.table(REPEAT.MX, file=paste0(out.dir, "/hg19repeats_repName.txt"),
            row.names=FALSE, col.names=TRUE, quote=FALSE, sep=";")

# rm(list=ls()); gc()