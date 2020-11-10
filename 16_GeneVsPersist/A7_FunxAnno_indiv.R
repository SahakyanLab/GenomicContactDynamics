################################################################################
# GO term and KEGG pathway enrichment for genes per Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/5_GeneVsPersist"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/3_AnnotationVsPersist"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
# Converter of HUGO gene symbols to ncbi-geneid for KEGG
hugoEntrezPath = paste0(data.dir, "/ucsc_tables/hsa_geneAnno/keyType_conversion/hg19anno_SYMBOLtoENTREZID_052020")
genelist.dir = paste0(wk.dir, "/out_anno_union")
out.dir = paste0(wk.dir, "/out_FunxAnno_indiv")
### OTHER SETTINGS #############################################################
# gcb + refseq
prefix = "min2Mb_ALL"
# Forebackground strings defining foreground and background, foreground;background
# to differentiate TLC from LC 
foreback.combi = c("cp_21;cp_HiC_all", "cp_21;cp_1") #"cp_21;cp_1") #"cp_21;cp_HiC_all" | "cp_21;genes" 'genes' for all genes
# Identifier of the foreground background combination. cp is treated as 
# id along with celltiss.v. No "cp" in 21 tissues/cellline so no conflicts.
# id.combi vector should correspond with foreback.combi.
id.combi = c("cp", "cp") #c("cp", "Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", 
             #  "SB", "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
plotOnly = FALSE
# If toPlot.v=NULL, faceted plot
toPlot.v = NULL
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
library(org.Hs.eg.db)
ORG = c(`GO_ALL`="org.Hs.eg.db", KEGG="hsa")
library(clusterProfiler)
library(ggplot2)
source(paste0(lib, "/GG_bgr.R"))
source(paste0(lib, "/funxAnno.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
genes.str <- list()
genes.str[["GO_ALL"]] <- readLines(con=paste0(genelist.dir, "/", prefix, "_name2"))
# Same order as genesEntrez
genes.str[["KEGG"]] <- readLines(con=paste0(genelist.dir, "/", prefix, "_entrezID"))
ind <- sapply(X=paste0("_", unlist(strsplit(x=foreback.combi, split=";")), "_end"), 
              FUN=grep, x=genes.str[["GO_ALL"]], fixed=TRUE)
ind <- unique( sort(c(ind, ind+1L), decreasing=FALSE) )
genes.str[["GO_ALL"]] <- genes.str[["GO_ALL"]][ind]
genes.str[["KEGG"]] <- genes.str[["KEGG"]][ind]; rm(ind)

foreback.combi.len <- length(foreback.combi)
if(length(id.combi)!=foreback.combi.len){
  stop("id.combi and foreback.combi not the same length.")
}
#---------------------------------------
# For converting ncbi-geneid to HUGO symbols for output
HE.conv <- read.delim(file=hugoEntrezPath, header=TRUE, row.names=NULL, 
                      stringsAsFactors=FALSE)
HE.conv <- HE.conv[!is.na(HE.conv$SYMBOL) & !is.na(HE.conv$ENTREZID),]
#HE.conv <- convertGeneKeyType(genes=NULL, convTablePath=hugoEntrezPath, drop.NA=TRUE)
HE.conv$ENTREZID <- as.character(HE.conv$ENTREZID)
#---------------------------------------
inputKey.v <- c(`GO_ALL`="SYMBOL", KEGG="ncbi-geneid")

for(i in 1:foreback.combi.len){
  forbgr.lab <- gsub(x=foreback.combi[i], pattern="\\;|\\_|end", 
                     replacement="")
  out.name <- paste0(prefix, "_", id.combi[i], "_", forbgr.lab)
  if(plotOnly==FALSE){
    
    # Get foreground and background genes from txt file
    # genes.lst[[1]] - foreground genes; genes.lst[[2]] - background genes
    funxAnnoOut <- list()
    for( appro in c("GO_ALL", "KEGG") ){
      genes.lst <- sapply(X=strsplit(x=foreback.combi[i], split=";")[[1]], simplify=FALSE, 
                          FUN=function(x){
                            # Extra "_" to differentiate LC from TLC and "_end" to differentiate
                            # _cp_1_end from _cp_10_end etc.
                            ind <- grep(pattern=paste0("_", x, "_end"), 
                                        x=genes.str[[appro]], fixed=TRUE)
                            if(length(ind)==1){
                              genes <- genes.str[[appro]][ind+1L]; rm(ind)
                            } else {
                              stop("Pattern not selective.")
                            }
                            # Also split with comma because entrezIds (when gene names are converted to
                            # entrezIds for KEGG) of one gene are separated by comma
                            genes <- strsplit(x=genes, 
                                              split=paste0("\\;|\\,"))[[1]]
                            unique(genes[!is.na(genes) & genes!="NA" & genes!="" & genes!=" "])
                          })
      funxAnnoOut[[appro]] <- funxAnno(input=genes.lst, org=ORG[[appro]], 
                                       inputKey=inputKey.v[appro], approach=appro, filePath=NULL)
      rm(genes.lst, appro); gc()
    }
    funxAnnoOut <- do.call("rbind", funxAnnoOut)
    rownames(funxAnnoOut) <- NULL
    
    # Convert ncbi geneids to HUGO
    KEGG.TF <- funxAnnoOut$ONTOLOGY=="KEGG"
    if( sum(KEGG.TF)>0 ){
      funxAnnoOut[KEGG.TF,"geneID"] <- sapply(X=funxAnnoOut[KEGG.TF, "geneID"], 
                                                   simplify=TRUE, FUN=function(id){
        id <- unique(strsplit(x=id, split="/", fixed=TRUE)[[1]])
        paste(x=HE.conv[HE.conv$ENTREZID%in%id,"SYMBOL"],
              collapse="/")
      })
    }
    rm(KEGG.TF); gc()
    
    write.csv(funxAnnoOut, file=paste0(out.dir, "/", out.name, ".csv"),
              row.names=FALSE, quote=FALSE)
    
  } else {
    funxAnnoOut <- read.csv(file=paste0(out.dir, "/", out.name, ".csv"),  
                            header=TRUE, stringsAsFactors=FALSE)
  }
  #---------------------------------------
  if(nrow(funxAnnoOut)>0){
    funxAnnoOut <- funxAnnoOut[order(funxAnnoOut$p.adjust),]
    funxAnnoOut$Description <- factor(funxAnnoOut$Description, 
                                      levels=unique(rev(funxAnnoOut$Description)))
    funxAnnoOut$Count <- as.numeric(as.character(funxAnnoOut$Count))
    
    if( is.null(toPlot.v) ){ toPlot.v <- unique(funxAnnoOut$ONTOLOGY) }
    ggplot( data=funxAnnoOut[funxAnnoOut$ONTOLOGY%in%toPlot.v,], 
            aes(x=-log10(p.adjust), y=Description) ) +
      geom_point(aes(colour=ONTOLOGY, size=Count)) + 
      geom_vline(xintercept=-log10(0.05), linetype="dashed", colour="red", size=0.7) +
      guides(colour="legend") +
      labs(title=paste0(out.name),
           colour=NULL, 
           size="Gene count",
           x=bquote(bold( "-log10("~"p-value"^"adj"~")" )), y=NULL
      ) +
      bgr5 + 
      facet_grid(~ONTOLOGY)
    
    ggsave(filename=paste0(out.dir, "/", out.name, "_", paste(sort(toPlot.v), collapse="_"), 
                           ".pdf"), units="in",
           width=15, height=20, limitsize=FALSE)
    
  } else {
    stop(paste0("No ", appro, " term enrichment result for", forbgr.lab, "."))
  }
  print(paste0(out.name, " done!"), quote=FALSE)
  rm(funxAnnoOut, forbgr.lab, out.name) ; gc()
} # foreback.combi.len end for loop

# rm(list=ls())
################################################################################

> genesSymbol[seq(from=1, to=340, by=2)]
[1] ">all_genes_end"                   
[2] ">all_genes_cp_HiC_all_end"        
[3] ">all_genes_cp_1_end"              
[4] ">all_genes_cp_2_end"              
[5] ">all_genes_cp_3_end"              
[6] ">all_genes_cp_4_end"              
[7] ">all_genes_cp_5_end"              
[8] ">all_genes_cp_6_end"              
[9] ">all_genes_cp_7_end"              
[10] ">all_genes_cp_8_end"              
[11] ">all_genes_cp_9_end"              
[12] ">all_genes_cp_10_end"             
[13] ">all_genes_cp_11_end"             
[14] ">all_genes_cp_12_end"             
[15] ">all_genes_cp_13_end"             
[16] ">all_genes_cp_14_end"             
[17] ">all_genes_cp_15_end"             
[18] ">all_genes_cp_16_end"             
[19] ">all_genes_cp_17_end"             
[20] ">all_genes_cp_18_end"             
[21] ">all_genes_cp_19_end"             
[22] ">all_genes_cp_20_end"             
[23] ">all_genes_cp_21_end"             
[24] ">all_genes_Co_count_HiC_all_end"  
[25] ">all_genes_Co_count_1_end"        
[26] ">all_genes_Co_count_2_end"        
[27] ">all_genes_Co_count_3_end"        
[28] ">all_genes_Co_count_4_end"        
[29] ">all_genes_Co_count_5_end"        
[30] ">all_genes_Co_count_grthan5_end"  
[31] ">all_genes_Hi_count_HiC_all_end"  
[32] ">all_genes_Hi_count_1_end"        
[33] ">all_genes_Hi_count_2_end"        
[34] ">all_genes_Hi_count_3_end"        
[35] ">all_genes_Hi_count_4_end"        
[36] ">all_genes_Hi_count_5_end"        
[37] ">all_genes_Hi_count_grthan5_end"  
[38] ">all_genes_Lu_count_HiC_all_end"  
[39] ">all_genes_Lu_count_1_end"        
[40] ">all_genes_Lu_count_2_end"        
[41] ">all_genes_Lu_count_3_end"        
[42] ">all_genes_Lu_count_4_end"        
[43] ">all_genes_Lu_count_5_end"        
[44] ">all_genes_Lu_count_grthan5_end"  
[45] ">all_genes_LV_count_HiC_all_end"  
[46] ">all_genes_LV_count_1_end"        
[47] ">all_genes_LV_count_2_end"        
[48] ">all_genes_LV_count_3_end"        
[49] ">all_genes_LV_count_4_end"        
[50] ">all_genes_LV_count_5_end"        
[51] ">all_genes_LV_count_grthan5_end"  
[52] ">all_genes_RV_count_HiC_all_end"  
[53] ">all_genes_RV_count_1_end"        
[54] ">all_genes_RV_count_2_end"        
[55] ">all_genes_RV_count_3_end"        
[56] ">all_genes_RV_count_4_end"        
[57] ">all_genes_RV_count_5_end"        
[58] ">all_genes_RV_count_grthan5_end"  
[59] ">all_genes_Ao_count_HiC_all_end"  
[60] ">all_genes_Ao_count_1_end"        
[61] ">all_genes_Ao_count_2_end"        
[62] ">all_genes_Ao_count_3_end"        
[63] ">all_genes_Ao_count_4_end"        
[64] ">all_genes_Ao_count_5_end"        
[65] ">all_genes_Ao_count_grthan5_end"  
[66] ">all_genes_PM_count_HiC_all_end"  
[67] ">all_genes_PM_count_1_end"        
[68] ">all_genes_PM_count_2_end"        
[69] ">all_genes_PM_count_3_end"        
[70] ">all_genes_PM_count_4_end"        
[71] ">all_genes_PM_count_5_end"        
[72] ">all_genes_PM_count_grthan5_end"  
[73] ">all_genes_Pa_count_HiC_all_end"  
[74] ">all_genes_Pa_count_1_end"        
[75] ">all_genes_Pa_count_2_end"        
[76] ">all_genes_Pa_count_3_end"        
[77] ">all_genes_Pa_count_4_end"        
[78] ">all_genes_Pa_count_5_end"        
[79] ">all_genes_Pa_count_grthan5_end"  
[80] ">all_genes_Sp_count_HiC_all_end"  
[81] ">all_genes_Sp_count_1_end"        
[82] ">all_genes_Sp_count_2_end"        
[83] ">all_genes_Sp_count_3_end"        
[84] ">all_genes_Sp_count_4_end"        
[85] ">all_genes_Sp_count_5_end"        
[86] ">all_genes_Sp_count_grthan5_end"  
[87] ">all_genes_Li_count_HiC_all_end"  
[88] ">all_genes_Li_count_1_end"        
[89] ">all_genes_Li_count_2_end"        
[90] ">all_genes_Li_count_3_end"        
[91] ">all_genes_Li_count_4_end"        
[92] ">all_genes_Li_count_5_end"        
[93] ">all_genes_Li_count_grthan5_end"  
[94] ">all_genes_SB_count_HiC_all_end"  
[95] ">all_genes_SB_count_1_end"        
[96] ">all_genes_SB_count_2_end"        
[97] ">all_genes_SB_count_3_end"        
[98] ">all_genes_SB_count_4_end"        
[99] ">all_genes_SB_count_5_end"        
[100] ">all_genes_SB_count_grthan5_end"  
[101] ">all_genes_AG_count_HiC_all_end"  
[102] ">all_genes_AG_count_1_end"        
[103] ">all_genes_AG_count_2_end"        
[104] ">all_genes_AG_count_3_end"        
[105] ">all_genes_AG_count_4_end"        
[106] ">all_genes_AG_count_5_end"        
[107] ">all_genes_AG_count_grthan5_end"  
[108] ">all_genes_Ov_count_HiC_all_end"  
[109] ">all_genes_Ov_count_1_end"        
[110] ">all_genes_Ov_count_2_end"        
[111] ">all_genes_Ov_count_3_end"        
[112] ">all_genes_Ov_count_4_end"        
[113] ">all_genes_Ov_count_5_end"        
[114] ">all_genes_Ov_count_grthan5_end"  
[115] ">all_genes_Bl_count_HiC_all_end"  
[116] ">all_genes_Bl_count_1_end"        
[117] ">all_genes_Bl_count_2_end"        
[118] ">all_genes_Bl_count_3_end"        
[119] ">all_genes_Bl_count_4_end"        
[120] ">all_genes_Bl_count_5_end"        
[121] ">all_genes_Bl_count_grthan5_end"  
[122] ">all_genes_MesC_count_HiC_all_end"
[123] ">all_genes_MesC_count_1_end"      
[124] ">all_genes_MesC_count_2_end"      
[125] ">all_genes_MesC_count_3_end"      
[126] ">all_genes_MesC_count_4_end"      
[127] ">all_genes_MesC_count_5_end"      
[128] ">all_genes_MesC_count_grthan5_end"
[129] ">all_genes_MSC_count_HiC_all_end" 
[130] ">all_genes_MSC_count_1_end"       
[131] ">all_genes_MSC_count_2_end"       
[132] ">all_genes_MSC_count_3_end"       
[133] ">all_genes_MSC_count_4_end"       
[134] ">all_genes_MSC_count_5_end"       
[135] ">all_genes_MSC_count_grthan5_end" 
[136] ">all_genes_NPC_count_HiC_all_end" 
[137] ">all_genes_NPC_count_1_end"       
[138] ">all_genes_NPC_count_2_end"       
[139] ">all_genes_NPC_count_3_end"       
[140] ">all_genes_NPC_count_4_end"       
[141] ">all_genes_NPC_count_5_end"       
[142] ">all_genes_NPC_count_grthan5_end" 
[143] ">all_genes_TLC_count_HiC_all_end" 
[144] ">all_genes_TLC_count_1_end"       
[145] ">all_genes_TLC_count_2_end"       
[146] ">all_genes_TLC_count_3_end"       
[147] ">all_genes_TLC_count_4_end"       
[148] ">all_genes_TLC_count_5_end"       
[149] ">all_genes_TLC_count_grthan5_end" 
[150] ">all_genes_ESC_count_HiC_all_end" 
[151] ">all_genes_ESC_count_1_end"       
[152] ">all_genes_ESC_count_2_end"       
[153] ">all_genes_ESC_count_3_end"       
[154] ">all_genes_ESC_count_4_end"       
[155] ">all_genes_ESC_count_5_end"       
[156] ">all_genes_ESC_count_grthan5_end" 
[157] ">all_genes_FC_count_HiC_all_end"  
[158] ">all_genes_FC_count_1_end"        
[159] ">all_genes_FC_count_2_end"        
[160] ">all_genes_FC_count_3_end"        
[161] ">all_genes_FC_count_4_end"        
[162] ">all_genes_FC_count_5_end"        
[163] ">all_genes_FC_count_grthan5_end"  
[164] ">all_genes_LC_count_HiC_all_end"  
[165] ">all_genes_LC_count_1_end"        
[166] ">all_genes_LC_count_2_end"        
[167] ">all_genes_LC_count_3_end"        
[168] ">all_genes_LC_count_4_end"        
[169] ">all_genes_LC_count_5_end"        
[170] ">all_genes_LC_count_grthan5_end"  


