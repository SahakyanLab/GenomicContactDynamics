################################################################################
# Fill upper triangle matrix using upper triangle df of contacts
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
# library(reshape2)
### FUNCTION ###################################################################
UpdfToUpTriMx <- function(df, mx.len='length/width of symmetrical matrix'){
  
  if( any(df$i>=df$j) ){
    stop("upperdfToMx(): Not upper triangle contacts.")
  }
  
  if( !is.numeric(df$value) ){
    stop("upperdfToMx(): Value not numeric.")
  }
  
  df <- df[order(df$j, df$i),]

  # Convert to matrix format
  MX <- matrix(data=NA, nrow=mx.len, ncol=mx.len)
  MX[ upper.tri(MX, diag=F) ] <- df$value

  return(MX)

}
################################################################################

# Filling upper matrix, originally I used Method 2 but shifted to shorter Method 1

# Starting upper triangle df from getContactDF(); not ordered
df.up <- expand.grid(i=1:5, j=1:5)
df.up <- df.base[df.up$i<df.up$j,]
df.up$value <- 1:10

# Method 1 (directly fills upper triangle)
df.1 <- df.up[order(df.up$j, df.up$i),]
MX <- matrix(data=NA, nrow=5, ncol=5)
MX[ upper.tri(MX, diag=F) ] <- df.1$value

# Method 2 (longer, fills lower triangle -> transpose to get upper triangle only)
df.2 <- df.up[order(df.up$i, df.up$j),]
MX[ lower.tri(MX, diag=F) ] <- df.2$value

df.1
df.2
MX
isSymmetric.matrix(MX)

> df.1
i j value
6  1 2     1
11 1 3     2
12 2 3     3
16 1 4     4
17 2 4     5
18 3 4     6
21 1 5     7
22 2 5     8
23 3 5     9
24 4 5    10
> df.2
i j value
6  1 2     1
11 1 3     2
16 1 4     4
21 1 5     7
12 2 3     3
17 2 4     5
22 2 5     8
18 3 4     6
23 3 5     9
24 4 5    10
> MX
[,1] [,2] [,3] [,4] [,5]
[1,]   NA    1    2    4    7
[2,]    1   NA    3    5    8
[3,]    2    3   NA    6    9
[4,]    4    5    6   NA   10
[5,]    7    8    9   10   NA
> isSymmetric.matrix(MX)
[1] TRUE


