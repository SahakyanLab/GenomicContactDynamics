lib = "/Users/ltamon/DPhil/lib"
source(paste0(lib, "/TrantoRextr/GEN_getKmers.R"))

library(Biostrings)
k.count <- getKmers(seq.string = DNAStringSet("TATATATATATAT"),  
                    k=7, method="Biostrings")
k.count[k.count > 0]

# When calculating k-mer counts
#1. Collapse sequence after masking removing

# same unmasked length, different k-mer counts  due to positioning of masked seq

"TATATATNNNAG"
"TATATATNAGNN"

"TATATATNNNAGCGG"
"TATATATNTATATNN"

# to make this equal?
"TATATATAGNNN"
"TATATATNAGNN"

#2. Normalise by unmasked length 

"TATATAT"
"TATATAT"

"TATATATATATAT"  
"TATATAT"       
ATATATA   TATATAT 
3/13      4/13 
          TATATAT 
          1/7
0.3956044

"TATATATATATATATA"
"TATATAT"
ATATATA TATATAT 
5/16    5/16 
        TATATAT 
        1/7
0.4821429

# For above case, measure sd of unmasked lengths because decreasing trend of sd
# can cause increasing complementarity

"TATATATATGCTAGTT"
"TATATAT"
ATATATA ATATATG ATATGCT ATGCTAG GCTAGTT TATATAT TATATGC TATGCTA TGCTAGT 
1       1       1       1       1       2       1       1       1 
                                        TATATAT
                                        1
0.5178571

#? Normalise by total k-mer count instead? 