################################################################################
RNAduplexR <- function(seq=seq, # vector of sequences to calculate data on, pairs
                       # should be next to each other, counting from the top
                       # Ex. c(seq1, seq2, seq3, seq4)
                       # >seq1
                       # >seq2
                       # .((.((.(((.((.&))))).))...)).   8,21  :   1,14  (-8.30)
                       # >seq3
                       # >seq4
                       # .((.((.(((.((.&))))).))...)).   8,21  :   1,14  (-8.30)
                       RNAduplex.CALL="/home/ltamon/prog/ViennaRNA-2.4.11/bin/RNAduplex",
                       # ID of temp file to avoid overwriting when in parallel
                       id.temp = "",
                       keep.tmp = FALSE,
                       NuclAc.type = "DNA",
                       filepath.par = "./dna_mathews.par"){
  
  # Check if each sequence has a pair 
  if( (length(seq)%%2)!=0 ){
    stop("A sequence is missing a pair.")
  }
  
  seq.len <- length(seq)
  
  temp.file <- paste0("temp3671236816_seq", id.temp, ".db")                   
  write(seq, file=temp.file)
  rm(seq); gc()
  
  # RNA parameters are described in
  # Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH. (2004)
  # Incorporating chemical modification constraints into a dynamic programming
  # algorithm for prediction of RNA secondary structure. Proc Natl Acad Sci USA 101(19):7287-92.

  # Ts will be treated as Us, but not converted into Us in the output,
  # though converted in the temp files, if preserved via keep.tmp.

  print(paste0("Folding ", NuclAc.type ,"s..."), quote=FALSE)
  #------------------------
  if(NuclAc.type == "RNA"){
    system(paste0(RNAduplex.CALL, 
                  " -d2 --noLP < temp3671236816_seq", id.temp, ".db > temp3671236816_rnaduplex", id.temp, ".txt")
           )
  }
  
  start.time <- Sys.time()
  
  #------------------------
  if(NuclAc.type == "DNA"){
    system(paste0(RNAduplex.CALL,
                  " --noconv -P ", filepath.par, " -d2 --noLP < temp3671236816_seq", id.temp, ".db > temp3671236816_rnaduplex", id.temp, ".txt")
           )
  }
  #------------------------
  
  end.time <- Sys.time()
  print(end.time-start.time)
  cat(end.time-start.time)

  OUTFILE <- readLines(paste0("temp3671236816_rnaduplex", id.temp, ".txt"))
  RNA.E <- RNA.STR.SCH <- RNA.FRQ <- RNA.PAIR.START <- RNA.PAIR.END <- rep(NA, seq.len/2)

  for( line in 1:(seq.len/2) ){
    
    print(line, quote=FALSE)
    
    splt <- unlist( strsplit(x=OUTFILE[line], split=" ") )
    
    RNA.E[line] <- as.numeric( strsplit(x=splt[11], split='\\(|\\)')[[1]][2] 
                               )
    RNA.STR.SCH[line] <- splt[1]
    
    splt4 <- strsplit(x=splt[4], split=",")[[1]]
    splt9 <- strsplit(x=splt[9], split=",")[[1]]
    
    RNA.PAIR.START <- paste(splt4[1], splt9[1], sep=";")
    RNA.PAIR.END <- paste(splt4[2], splt9[2], sep=";")
  
  }

  RESULT <- NULL
  RESULT$RNA.E          <- RNA.E
  RESULT$RNA.STR.SCH    <- RNA.STR.SCH
  RESULT$RNA.PAIR.START <- RNA.PAIR.START
  RESULT$RNA.PAIR.END   <- RNA.PAIR.END
 
  if(keep.tmp==FALSE){
    file.remove(c(paste0("temp3671236816_seq", id.temp, ".db"), 
                  paste0("temp3671236816_rnaduplex", id.temp, ".txt")
                ))
  }

  return( data.frame(RESULT, stringsAsFactors=FALSE) )

}
################################################################################
