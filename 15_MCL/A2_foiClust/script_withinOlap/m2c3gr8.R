################################################################################
# Mark nodes/vertices based on overlap with a feature then compare Cs and Cp
# results.
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelCluster" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/21_MCL"
    data.dir = "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/21_MCL"
    data.dir = "/t1-data/user/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
mcl.dir = paste0(wk.dir, "/out_mcl")
out.dir = paste0(wk.dir, "/out_foiClust")
foi.dir = paste0(wk.dir, "/foi")
#foi.dir = paste0(data.dir, "/funx_data_fixCoordSys/masterpool_hg19_convTo1based/reduced_b2")
foifile = paste0(wk.dir, "/foifile/foifile_FC")
hub.dir = paste0(wk.dir, "/out_hubfile")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = "chr3" #paste("chr", c(1:22, "X"), sep="")
bin.len = 40000
topCP = 21; cp.v = 1:21
ct = "FC"; ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
                    "AG", "Ov", "Bl", "MesC","MSC", "NPC", "TLC", "ESC", "FC", "LC")
gap.bin = 0
e.val = 2 # MCL - Expansion
gran.val = 8 #c(2,3,4,5,10,20,30,40,50,60,80) # MCL - Inflation

hub.id = "FC_topCP3_gapBin50"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(sna)
library(GenomicRanges)
library(yarrr)
source(paste0(lib, "/TrantoRextr/GEN_WhichOverlap.R"))
source(paste0(lib, "/finaliseFOI.R"))
source(paste0(wk.dir, "/lib/makeClustgplot.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
for(chr in chr.v){
  
chr.id <- paste0(gcb, "_", chr)
print(paste0(chr.id, "..."), quote=FALSE)

ct.id <- NULL
if( !is.null(ct) ){ ct.id <- ct; ct.id <- paste0(ct.id[!is.na(ct.id)], "_") }
out.id <- paste0(ct.id, "topCP", topCP, "_gapBin", gap.bin)  
mcl.id <- paste0("gran", gran.val, "_expan", e.val)
out.name <- paste0(chr.id, "_", out.id, "_", mcl.id)
Cs.id <- paste0(out.id, "_weightCs"); Cs.dir <- paste0(mcl.dir, "/", Cs.id, "_", mcl.id)
Cp.id <- paste0(out.id, "_weightCp"); Cp.dir <- paste0(mcl.dir, "/", Cp.id, "_", mcl.id)

# List of features
foi.v <- finaliseFOI(foi.dir=foi.dir, foifile=foifile)

M.lst <- sapply(X=c("Cs", "Cp"), simplify=FALSE, FUN=function(weight){
  eval(parse(text=paste0(
    'load(file=paste0(', weight, '.dir, "/", chr.id, "_",', weight, 
    '.id, "_", mcl.id, "_mcl.RData"))'
  )))
  return(MCLOUT)
})

PMX.lst <- sapply(X=c("Cs", "Cp"), simplify=FALSE, FUN=function(weight){
  eval(parse(text=paste0(
    'load(file=paste0(', weight, '.dir, "/", chr.id, "_",', weight, 
    '.id, "_", mcl.id, "_posmx.RData"))'
  )))
  return(POS.MX)
})

# Get list of nodes
node.v <- as.numeric(colnames(M.lst$Cs$equilibrium.state))
node.v.len <- length(node.v)
node.end <- node.v*bin.len
node.start <- node.end-bin.len+1
node.v <- as.character(node.v)
#-------------------------------------------------------------------------------
pdf(file=paste0(out.dir, "/", out.name, "_mclCsVsCp.pdf"), height=20, width=40)
par( mfcol=c(2,length(foi.v)+2) )
#---------------------------------------Colour by Cs cluster
# Assign colour to nodes based on Cs clusters
clust.col <- as.character(M.lst$Cs$cluster)
temp.clust <- setdiff(unique(clust.col),"0")
temp.col <- c("black", colorRampPalette(yarrr::piratepal("basel"))(length(temp.clust)))
names(temp.col) <- c("0", temp.clust)
clust.col <- adjustcolor(temp.col[clust.col], alpha.f=0.5) 
rm(temp.clust, temp.col)

for(weight in c("Cs", "Cp")){
  # Plot
  makeClustgplot(MCLeq.mx=M.lst[[weight]][["equilibrium.state"]], 
                 coord.mx=PMX.lst[[weight]], label=node.v, v.col=clust.col,
                 title=paste0(out.name, "_weight", weight, "_CsCluster"))
}
rm(clust.col)
print("Colouring by Cs cluster done!", quote=FALSE)
#---------------------------------------Colour by feature
for(foi in foi.v){
  
  foi.bed <- read.delim(file=paste0(foi.dir, "/", foi), stringsAsFactors=FALSE, 
                        header=FALSE)[,1:3]
  foi <- tail(x=strsplit(x=foi, split="\\/")[[1]], n=1)
  foi <- gsub(x=foi, pattern="ct\\_|foi\\_|desc\\_|\\.bed", replacement="")
  
  # Assign colour to nodes based on foi
  chr.TF <- foi.bed[,1]==chr
  hits.mx <- cbind(query=NA, subject=NA)
  if( sum(chr.TF)!=0 ){
    hits.mx <- WhichOverlap(start.query=node.start, 
                            end.query=node.end, 
                            space.query=chr,
                            start.subject=foi.bed[chr.TF,2],
                            end.subject=foi.bed[chr.TF,3],
                            space.subject=chr, 
                            maxgap=-1L, minoverlap=1L, type="within")
  } else{
    print(paste0(chr, ": No ", foi, "."), quote=FALSE)
  }
  hits.node <- node.v[ hits.mx[,"query"] ]; rm(hits.mx)
  
  node.col <- rep("black", times=node.v.len)
  names(node.col) <- node.v
  
  dup.TF <- duplicated(hits.node)
  if( any(dup.TF) ){
    print("Promiscuous nodes.", quote=FALSE)
    node.col[node.v%in%hits.node[dup.TF]] <- "green"
  }; rm(dup.TF)
  
  node.col[node.v%in%hits.node & node.col!="green"] <- "red"
  node.col <- adjustcolor(node.col, alpha.f=0.5)
  rm(chr.TF, hits.node); gc()
  
  for(weight in c("Cs", "Cp")){
    # Plot
    makeClustgplot(MCLeq.mx=M.lst[[weight]][["equilibrium.state"]], 
                   coord.mx=PMX.lst[[weight]], label=node.v, v.col=node.col,
                   title=paste0(out.name, "_weight", weight, "_", foi))
  }
  
  print(paste0(foi, " done!"), quote=FALSE)
  rm(foi.bed, foi, node.col); gc()
  
} # foi.v for loop end
print("Colouring by feature done!", quote=FALSE)
#---------------------------------------Colour by hub
hub.v <- list.files(path=hub.dir, pattern=paste0(chr.id, "_", hub.id, "_hub"))
temp.col <- colorRampPalette(yarrr::piratepal("basel"))(length(hub.v))
names(temp.col) <- hub.v
hub.col <- rep("black", times=node.v.len)
for(hub in hub.v){
  hnode.v <- as.character(read.csv(file=paste0(hub.dir, "/", hub), header=TRUE,
                                   stringsAsFactors=FALSE)[,"bin"])
  nonh.TF <- hub.col=="black"
  h.TF <- node.v%in%hnode.v
  hub.col[h.TF & nonh.TF] <- temp.col[hub]
  # Promiscous node (node shared by >1 hubs)
  hub.col[h.TF & !nonh.TF] <- "green"
  rm(nonh.TF, h.TF)
}
v.cex <- rep(1, times=node.v.len); v.cex[hub.col!="black"] <- 2
hub.col <- adjustcolor(hub.col, alpha.f=0.5)
for(weight in c("Cs", "Cp")){
  # Plot
  makeClustgplot(MCLeq.mx=M.lst[[weight]][["equilibrium.state"]], 
                 coord.mx=PMX.lst[[weight]], label=node.v, v.col=hub.col, v.cex=v.cex,
                 title=paste0(out.name, "_weight", weight, "_hub_", hub.id))
}
rm(hub.col)
print("Colouring by hub done!", quote=FALSE)
#---------------------------------------
dev.off()
print(paste0(chr, " done!"), quote=FALSE)

}

# rm(list=ls()); gc()






