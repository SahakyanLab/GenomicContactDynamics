################################################################################
# Test edlib alignment modes
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    wk.dir = "/Users/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir =  "/Users/ltamon/Database"
  } else if(whorunsit == "LiezelCluster"){
    lib = "/t1-data/user/ltamon/DPhil/lib"
    wk.dir = "/t1-data/user/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir =  "/t1-data/user/ltamon/Database"
  } else if(whorunsit == "LiezelLinuxDesk"){
    lib = "/home/ltamon/DPhil/lib"
    wk.dir = "/home/ltamon/DPhil/GenomicContactDynamics/11_Constraints"
    data.dir =  "/home/ltamon/Database"
  } else {
    print("The supplied <whorunsit> option is not created in the script.", quote=FALSE)
  }
}
### OTHER SETTINGS #############################################################
type = "align-HW"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(Rcpp)
if(type=="align-NW"){
  
  print("align-NW", quote=FALSE)
  sourceCpp(file=paste0(wk.dir, "/lib/edlibNW.cpp"))
  
} else if(type=="align-HW"){
  
  print("align-HW", quote=FALSE)
  sourceCpp(file=paste0(wk.dir, "/lib/edlibHW.cpp"))
  
} else {
  
  stop("Invalid input for type.")
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# EDLIB_MODE_HW - deleting ends of 2nd sequence (target) is free but not of 1st
# sequence (query); substitution penalised with 1

# This will result to a biased treatment of sequences depending which is the target/
# query (given that the lengths are different). Therefore, two permutations 
# should be tested. Ex.
BBBC-
---CA
11100=3

---CA
BBBC-
00001=1

#------------------------------------------------------------------------------

# EDLIB_MODE_NW - gap/substitution penalised with 1

# EDLIB_MODE_NW does not only align sequence from end to end and then
# calculate edit distance. It matches the two sequences looking for
# common patches, introducing gaps along the way. Then both the gap
# and substitution is penalised 1. Edit distance will be the same regardless
# of the query/target assignment

#---------------------------------------
# 6
edlibNW(query="TELEPHONE",
            target="LEP", 
            queryLength=9,
            targetLength=3)
        
TELEPHONE
--LEP----
110001111=6  

#---------------------------------------
# 6
edlibNW(query="TELEPHONE",
        target="ALEP", 
        queryLength=9,
        targetLength=4)
TELEPHONE
-ALEP----
110001111=6  
#---------------------------------------
# 9
edlibNW(query="TELEPHONE",
        target="AAAAAAAAA", 
        queryLength=9,
        targetLength=9)
TELEPHONE
AAAAAAAAA
111111111=9 
#---------------------------------------
# 9
edlibNW(query="TELEPHONE",
        target="", 
        queryLength=9,
        targetLength=0)
TELEPHONE
---------
111111111=9 
#---------------------------------------
# 8
edlibNW(query="TELEPHONE",
        target="AAAALEPAA", 
        queryLength=9,
        targetLength=9)
--TELEPHONE
AAAALEPAA
11110001111=8
#---------------------------------------
# 9
edlibNW(query="TELEPHONE",
        target="AAAALFPAA", 
        queryLength=9,
        targetLength=9)
--TELEPHONE
AAAALFPAA
11110101111=9
#---------------------------------------
# 6
edlibNW(query="TELEPHONE",
        target="ALEAAONPA", 
        queryLength=9,
        targetLength=9)
TELEPHONE
 ALEAAONPA
1100110011=6
#---------------------------------------

# EDLIB_MODE_NW over EDLIB_MODE_HW
# 1. Scoring results to greater difference in complementarity between contacts
# based on similar sequences they contain

################################################################################

# rm(list=ls())