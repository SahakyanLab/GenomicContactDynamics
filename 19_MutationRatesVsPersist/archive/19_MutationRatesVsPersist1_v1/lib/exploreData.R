################################################################################
# Explore COSMIC non-coding variant dataset.
### FUNCTION ###################################################################
exploreData <- function(df){
  
  x <- list()
  
  # Total mutations
  x$totmut <- length(df[,1])
  
  x$columns <- colnames(df)
  
  # Check for NAs 
  tmp <- list()
  for( i in 1:ncol(df) ){
    
    colnme <- colnames(df)[i]
    NA.TF <- is.na(df[,colnme])
    if( any(NA.TF) ){
      tmp[[colnme]] <- colnme
    }
    #df[NA.TF,colnme] <- "Missing"
    rm(colnme)
    
  }
  x$withNAs <- names(tmp)
  rm(tmp)
  
  # ID_SAMPLE
  x$ID_SAMPLE <- table(ncv.df$ID_SAMPLE, useNA="always")
  
  # Sample name
  x$`Sample name`<- table(ncv.df$`Sample name`, useNA="always")
  
  # Primary site
  x$`Primary site` <- table(df$`Primary site`, useNA="always")
  
  # Primary histology
  x$`Primary histology` <- table(df$`Primary histology`, useNA="always")
  
  # Zygosity
  x$zygosity <- table(df$zygosity, useNA="always")
  
  # Genome reference
  x$GRCh <- table(df$GRCh, useNA="always")
  
  # Genomic position
  WT_SEQ.len <- nchar(df$WT_SEQ)
  MUT_SEQ.len <- nchar(df$MUT_SEQ)
  #SEQ.len <- pmax(ref.len, alt.len)

  position <- strsplit(x=df$`genome position`, split=":|-")
  
  pos.class <- unlist(lapply( X=1:length(position), FUN=function(i){
    
    x <- position[[i]]
    
    if( any( is.na(x) | x=="NA" ) ){
      return(NA)
    } else if( length(x==3) ){
      
      x <- as.numeric(x)
      if( (x[3]-x[2]+1L)==WT_SEQ.len[i] ){
        return("onebased")
      } else if( (x[3]-x[2])==WT_SEQ.len[i] ){
        return("zerobased")
      } else {
        return("weird")
      }
      
    } else {
      stop("Wrong genome position.")
    }
    
  }))
  rm(i, SEQ.len)
  
  x$position <- table(as.vector(pos.class), useNA="always")
  
  # Chromosome
  x$chr <- unlist(lapply(X=position, FUN=function(x) x[1])) 
  x$chr <- table(x$chr, useNA="always")
 
  rm(position)
  
  # Mutation somatic status
  x$`Mutation somatic status` <- table(df$`Mutation somatic status`, useNA="always")
  
  # WT_SEQ and MUT_SEQ
  x$WT_SEQ <- table(df$WT_SEQ, useNA="always")
  x$MUT_SEQ <- table(df$MUT_SEQ, useNA="always")
  
  df$WT_SEQ <- toupper(df$WT_SEQ)
  df$MUT_SEQ <- toupper(df$MUT_SEQ)
  base.v <- c("A", "C", "G", "T")
  #base.v <- c(base.v, tolower(base.v))
  x$validSingleBaseMut <- sum(nchar(df$WT_SEQ)==1 & MUT_SEQ.len==1 & 
                              df$WT_SEQ%in%base.v & df$MUT_SEQ%in%base.v & 
                              df$WT_SEQ!=df$MUT_SEQ &
                              !is.na(df$`genome position`))
  rm(base.v)
  
  # SNP
  x$SNP <- table(df$SNP, useNA="always")
  
  # FATHMM_MKL_NON_CODING_SCORE & FATHMM_MKL_CODING_SCORE
  for( y in c("FATHMM_MKL_NON_CODING_SCORE", "FATHMM_MKL_CODING_SCORE") ){
    
    tmp <- df[[y]]
    NA.TF <- is.na(tmp)
    tmp <- tmp[!NA.TF]
    x[[y]] <- c(Neutral=sum(tmp<=0.5),
                Pathogenic=sum(tmp>=0.7), 
                Mid=sum(tmp>0.5 & tmp<0.7),
                Missing=sum(NA.TF))
    
    rm(tmp, y, NA.TF)
    
  }
  
  # FATHMM_MKL_NON_CODING_GROUPS & FATHMM_MKL_CODING_GROUPS
  for( y in c("FATHMM_MKL_NON_CODING_GROUPS", "FATHMM_MKL_CODING_GROUPS") ){
    
    x[[y]] <- table(df[[y]], useNA="always")
    rm(y)
    
  }
  
  # Whole Genome Reseq
  x$Whole_Genome_Reseq <- table(df$Whole_Genome_Reseq, useNA="always")
  
  # Whole_Exome
  x$Whole_Exome <- table(df$Whole_Exome, useNA="always")
  
  # Check
  nme.TF <- !names(x)%in%c("totmut", "columns", "withNAs")
  for( nme in names(x)[nme.TF] ){
    
    if(nme!="validSingleBaseMut"){
      if( sum(x[[nme]])!=x$totmut ){
        stop(paste0(nme, ": Sum not equal to total mutations."))
      }
    }
    
    x[[nme]] <- x[[nme]]/x$totmut
    
  }
  
  return(x)

}

# rm(list=ls()); gc()

