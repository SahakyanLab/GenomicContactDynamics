################################################################################
# Network representation of persistent substructures formed by highly persistent
# contacts. Added feature of marking areas overlapping with a feature provided
# as a bed file. 
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
data.dir = paste0(home.dir, "/Database")
wk.dir = paste0(home.dir, "/SahakyanLab/GenomicContactDynamics/23_Hub_Gene_Expression")
persist.dir = paste0(wk.dir, "/out_basePersist")
out.dir = paste0(wk.dir, "/out_network_test")
centrobed.file = paste0(wk.dir, "/txTable/ct_hg19_foi_centromoreonly_desc_DNA")
chrLenfile = paste0(data.dir, "/genome_info/Hsa_GRCh37_73_chr_info.txt")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chr.v = "chr21" #paste0("chr", c(1:22, "X"), sep="")
bin.len = 40000
# topCP=3 would means the top 3 values of cp.v; topCP=-3 means the bottom 3 values
# of cp.v; cp.v shouldn't have duplicates as better to arrange in increasing order
topCP = 3; cp.v = 1:21
ct = NULL; ct.v = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
                    "AG", "Ov", "Bl", "MesC","MSC", "NPC", "TLC", "ESC", "FC", "LC")
# Gap in terms of % of chr length ("Perc") or bin ("Bin")
gap.type = "Bin" 
gap.val = 200
# Number of bins in between non-contact bins (filler nodes)
edge.res = 50
# Multiplier to length of edges to make the difference in lengths more obvious;
# default=2L
edgelen.mult = 2L
#-------------------
marking = FALSE
# Bed file format (chr, start, end). Ranges will not be reduced to make them
# non-overlapping. Ideally, each range should corrrespond to a unique feature/marker.
# Internally, row numbers are added as uniqueID of the feature (before removing
# ranges with NA start and end coordinates to make sure that the row numbers
# correspond to the original bed file provided.)
#markfilePath = paste0(wk.dir, "/sample_markfile/sample_markfile_continuous")
markfilePath = paste0(wk.dir, "/A2_app_cge/sample_feature_bed/sample_markfile_discrete")
bed.CS = "0-based" # 1-based or 0-based
mark.id = "discrete_test" #"FC_Schmitt_TADboundary_trash"
header.bed = FALSE
colouring.style = "discrete"
# Add extra bins from markfile?
addNodes = FALSE
#-------------------
plotOnly = FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(RColorBrewer)
library(visNetwork)
library(compiler)
library(GenomicRanges)
library(dplyr)
library(viridis)
source(paste0(wk.dir, "/lib/GEN_WhichOverlap.R"))
source(paste0(wk.dir, "/lib/markCentromereOnNetworkData.R"))
source(paste0(wk.dir, "/lib/makeNetworkData.R"))
source(paste0(wk.dir, "/lib/modifyNetworkData.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Chromosome length file
chrLen.df <- read.table(file=chrLenfile, as.is=FALSE, header=TRUE,
                        colClasses=c("character", "integer", "integer"))
#---------------------------------------
cp.v <- sort(cp.v, decreasing=FALSE)
if(topCP<0){ cp.v <- rev(cp.v) }
cp <- sort(tail(cp.v, n=abs(topCP)), decreasing=FALSE)
#---------------------------------------
if(marking){
  markfile <- read.delim(file=markfilePath, header=header.bed, stringsAsFactors=FALSE)
}
#---------------------------------------
ct.id <- NULL
if( !is.null(ct) ){ ct.id <- ct; ct.id <- paste0(ct.id[!is.na(ct.id)], "_") }
#---------------------------------------
for(chr in chr.v){
  
  out.name <- paste0(gcb, "_", chr, "_", ct.id, "topCP", topCP, "_gap", gap.type, 
                     gap.val, "_edgeRes", edge.res, "_edgelenmult", edgelen.mult)
  if(marking){ out.name <- paste0(out.name, "_", mark.id) }
  
  if(plotOnly==FALSE){
    
    print("Regenerating plot data...", quote=FALSE)
   
    # Load PERSIST.MX; filter contacts based on Cp, ct and gap percentage threshold
    #load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, ".RData"))
    load(file=paste0(persist.dir, "/", chr, "_Persist_", gcb, "_topCP4.RData"))
    chr.len <- chrLen.df[chrLen.df$chromosome==chr,"length.bp"]
    
    gap.bin.v <- (PERSIST.MX$hits$j-PERSIST.MX$hits$i)-1
    if(gap.type=="Perc"){
      gap.TF <- (100*(gap.bin.v*bin.len/chr.len)) >= gap.val
    } else if(gap.type=="Bin"){
      gap.TF <- gap.bin.v >= gap.val
    } else {
      stop("Invalid gap.type (Bin/Perc only)") 
    }
    
    if( is.null(ct) ){
      incl.TF <- PERSIST.MX$ntis%in%cp & gap.TF
    } else if(ct%in%ct.v){
      incl.TF <- PERSIST.MX$ntis%in%cp & gap.TF & PERSIST.MX$hits[[ct]]>0
    } else {
      stop("Invalid ct argument. It should be NULL or an element of ct.v")
    }
    rm(gap.TF, gap.bin.v)
    if(sum(incl.TF)==0){
      print(paste0(chr, ": No contacts to be shown."), quote=FALSE)
    }
    #---------------------------------------Last bin of chr
    bin.last <- ceiling(chr.len/bin.len)
    if(bin.last!=chrLen.df[chrLen.df$chromosome==chr,"bins.40kb"]){
      stop("Incorrect last bin number.")
    }
    rm(chr.len); gc()
    #---------------------------------------
    mark.df <- NULL
    if(marking){
      mark.TF <- markfile[,1]==chr
      if(sum(mark.TF)>0){
        mark.df <- markfile[mark.TF,]
      }
    }
    NETWRK <- makeNetworkData(contact.mx=PERSIST.MX$hits[incl.TF,c("i","j")],
                              chr=chr,
                              cp.v=PERSIST.MX$ntis[incl.TF],
                              bin.last=bin.last,
                              edge.res=edge.res,
                              edgelen.mult=edgelen.mult,
                              mark.df=mark.df, 
                              bed.CS=bed.CS,
                              header.bed=header.bed,
                              olap.col=olap.col,
                              addNodes=addNodes,
                              centrobed.file=NULL)
    rm(PERSIST.MX); gc()
    save(NETWRK, file=paste0(out.dir, "/", out.name, ".RData"))
  } else {
    file.nme <- paste0(out.dir, "/", out.name, ".RData")
    if( file.exists(file.nme) ){
      load(file=file.nme)
    } else {
      print(paste0(chr, " NETWRK data does not exist!"), quote=FALSE)
      next
    }
  }
  #---------------------------------------Plot
  graph <- visNetwork(NETWRK$nodes, NETWRK$edges, 
                      main=list(text=paste0(out.name, "_Nij=", sum(!is.na(NETWRK$edges$label))),
                                style= "font-family:Georgia, Times New Roman, Times, serif;font-weight:bold;font-size:5px;text-align:center;"),
                      width=1600, height=1000
                      ) %>%
    # Values are the defaults
    visPhysics(solver="barnesHut", barnesHut=list(springConstant=0.04), 
               stabilization=list(enabled=TRUE), timestep=0.5) %>%
    visExport(type="pdf", name=out.name, label=paste0("Export as pdf"), float="right") %>%
    visSave(file=paste0(out.dir, "/", out.name, ".html"), selfcontained=TRUE, 
            background="white")
  rm(graph, NETWRK); gc()
  print(paste0(chr, " done!"), quote=FALSE)
}

# rm(list=ls()); gc()


