# "primary_cohort"
# "dixon_2012"
# "fraser_2016"
# "seitan_2013"
# "shen_2012"
# "sofueva_2013"
# "zuin_2014" 

#---------------------------------------------- find out the available file sets
folder = "zuin_2014"

cell.type <- sapply(dir(folder), FUN=function(i){
          strsplit(strsplit(i,".chr", fixed=T)[[1]][1], ".nor", fixed=T)[[1]][1]
                    })

chr <- sapply(dir(folder), FUN=function(i){
  strsplit(strsplit(i,".chr", fixed=T)[[1]][2], ".qq.mat", fixed=T)[[1]][1]
})

paste(cell.type,chr)
unique(cell.type)
#-------------------------------------------------------------------------------

#---------------------------------------- find out non-standard chromosomal sets
data.str <- read.table("data_structure.txt", as.is=TRUE, header=TRUE)
for(i in 1:dim(data.str)[1]){
  files <- grep( paste0(data.str[i,"cell_type"],".nor."),
                  dir(data.str[i,"source"]), value=TRUE )
  chrs <- sapply(files, FUN=function(fls){
    strsplit(strsplit(fls, ".chr", fixed=T)[[1]][2],".qq.mat",fixed=T)[[1]][1]
    }, USE.NAMES=F)
  
  mtch <- match( c(as.character(1:22),"X"), chrs)
  if(sum(is.na(mtch))==0 & length(c(as.character(1:22),"X"))==length(chrs)){
    print(paste(i,"human", sep=" "), quote=F)
  } else {
    print(paste(i, paste(chrs,collapse=" "),sep=" "), quote=F)
  }
}
#-------------------------------------------------------------------------------

# source cell_type merge_source_if_repl
# primary_cohort AD -
# primary_cohort AO -
# primary_cohort BL -
# primary_cohort CO -
# primary_cohort GM12878 -
# primary_cohort h1.merge h1.rep1;h1.rep2
# primary_cohort HC -
# primary_cohort IMR90 -
# primary_cohort LG.merge LG1;LG2
# primary_cohort LI -
# primary_cohort LV -
# primary_cohort mes.merge mes.rep1;mes.rep2
# primary_cohort msc.merge msc.rep1;msc.rep2
# primary_cohort npc.merge npc.rep1;npc.rep2
# primary_cohort OV -
# primary_cohort PA.merge PA2;PA3
# primary_cohort PO.merge PO1;PO3
# primary_cohort RV -
# primary_cohort SB2 -
# primary_cohort SX.merge SX1;SX3
# primary_cohort tro.merge tro.rep1;tro.rep2
# dixon_2012 mES.merge mES.rep1.HindIII;mES.rep2.HindIII
# dixon_2012 mES.NcoI -
# fraser_2016 HindIII_mESC -
# fraser_2016 HindIII_Neuron -
# fraser_2016 HindIII_NPC -
# fraser_2016 NcoI_mESC -
# fraser_2016 NcoI_Neuron -
# fraser_2016 NcoI_NPC -
# seitan_2013 Tcell.Rad21.KO.merge Tcell.Rad21.KO.rep1;Tcell.Rad21.KO.rep3
# seitan_2013 Tcell.Rad21.WT.merge Tcell.Rad21.WT.rep1;Tcell.Rad21.WT.rep3
# shen_2012 mCO.merge mCO.rep1.HindIII;mCO.rep2.HindIII
# sofueva_2013 mAST-Rad21-deleted.merge mAST-Rad21-deleted.rep1;mAST-Rad21-deleted.rep2
# sofueva_2013 mAST-Rad21-floxed.merge mAST-Rad21-floxed.rep1;mAST-Rad21-floxed.rep2
# sofueva_2013 mNSC-Rad21-deleted.merge mNSC-Rad21-deleted.rep1;mNSC-Rad21-deleted.rep2
# sofueva_2013 mNSC-Rad21-floxed.merge mNSC-Rad21-floxed.rep1;mNSC-Rad21-floxed.rep2
# sofueva_2013 mAST-Rad21-WT.AdCre.rep1 -
# sofueva_2013 mAST-Rad21-WT.rep1 -
# sofueva_2013 mNSC-Rad21-floxed-G1.rep1 -
# sofueva_2013 mNSC-Rad21-WT.OHT.rep1 -
# sofueva_2013 mNSC-Rad21-WT.rep1 -
# zuin_2014 HRV.merge HRV.rep1;HRV.rep2
# zuin_2014 siCONTROL.merge siCONTROL.rep1;siCONTROL.rep2 
# zuin_2014 siCTCF.merge siCTCF.rep1;siCTCF.rep2
# zuin_2014 TEV.merge TEV.rep1;TEV.rep2
