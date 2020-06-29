################################################################################
# Function to perform Markov Clustering of SPARSE graphs (<< number of edges 
# compared to highest possible number of edges) 
################################################################################
# LIBRARIES & DEPENDANCES * LIBRARIES & DEPENDANCIES * LIBRARIES & DEPENDANCES *
################################################################################
#library(expm)
### FUNCTION ###################################################################
MyMCL <- function(adjMX, e=2, gran=2, max.iter=1000, 
                  # Don't add self-loops if FALSE
                  addLoops=1,
                  # Allow Vertices to form own Cluster? If
                  allowVC=NULL) {
  if( is.null(addLoops) ){ stop("Add self loops? Value or FALSE.") }
  if( is.null(allowVC) ){ stop("Allow vertices to form own cluster? TRUE/FALSE.") }
  if( !is.null(addLoops) ){
    diag(adjMX) <- addLoops
  } else if(addLoops==FALSE){
    print("Not adding self loops", quote=FALSE)
  }   
  dimn.input <- dimnames(adjMX)
  #-------------------Apply MCL algorithm; expansion-inflation
  # Initial scaling
  adjMX <- apply( X=adjMX, MARGIN=2, FUN=function(col){col/sum(col)} )
  iter <- 1
  repeat{
    #adjMX1 <- adjMX%*%adjMX # Expansion
    adjMX1 <- adjMX%^%e   # Expansion; operator from expm library
    adjMX1 <- adjMX1^gran # Inflation 
    adjMX1 <- apply( X=adjMX1, MARGIN=2, FUN=function(col){col/sum(col)} ) # Scale
    if( identical(adjMX,adjMX1) ) {
      rm(adjMX1)
      print("MCL iteration")
      print(paste0("CONVERGED AFTER ", iter, " iterations!")); break 
    } else if(iter==max.iter){
      stop("DID NOT CONVERGED. INCREASE MAX ITERATION.")
    } else {
      print("MCL iteration"); adjMX <- adjMX1
      iter <- iter+1
    }
  }
  #-------------------Identify clusters (Fom "MCL" R lib (2015) with modifications)
  node.num <- ncol(adjMX) 
  clust <- matrix(adjMX[rowSums(adjMX)!=0,], ncol=node.num)
  attr.num <- nrow(clust) # Number of attractors (nodes with non-zero rowSum)
  if(attr.num==0){ stop("No attractors!") }
  for(i in 1:attr.num){
    for(j in 1:node.num) {
      if((clust[i,j] < 1) & (clust[i,j] > 0)){
        clust[,j] <- 0
        clust[i,j] <- 1
      }
    }
  }
  for(i in 1:attr.num){ clust[i,][clust[i,]!=0] <- i }
  clust.grp <- colSums(clust)
  rm(attr.num)
  #-------------------Allow vertices to form own cluster?
  if(!allowVC){
    dup <- duplicated(clust.grp) + duplicated(clust.grp, fromLast=T)
    clust.grp[dup==0] <- 0; rm(dup)
  }
  #-------------------
  if( is.null(dimn.input) | is.null(dimn.input[[1]]) | is.null(dimn.input[[2]]) ){
    dimnames(adjMX) <- list(1:node.num, 1:node.num); rm(node.num)
  } else {
    dimnames(adjMX) <- dimn.input
  }
  OUTPUT <- list(K=sum(unique(clust.grp)!=0),
                 converged.iter=iter,
                 cluster=unname(clust.grp),
                 equilibrium.state=adjMX)
  return(OUTPUT)
}

#Sample: 5 nodes, 3 edges
#Cluster A: 1-2-5  
#Cluster B: 3-4

# In the MCL converged output below, cluster A represented
# as 1s at the central node, 2, while cluster B (pair)
# is represented as 0.5 at 3 and 4 (since no central node). 
#> adjMX1
#      [,1] [,2] [,3] [,4] [,5]
#[1,]    0    0  0.0  0.0    0
#[2,]    1    1  0.0  0.0    1
#[3,]    0    0  0.5  0.5    0
#[4,]    0    0  0.5  0.5    0
#[5,]    0    0  0.0  0.0    0

# To deduce cluster, get clust matrix, containing rows of adjMX1
# with rowSum=0
#> clust
#      [,1] [,2] [,3] [,4] [,5]
#[1,]    1    1  0.0  0.0    1
#[2,]    0    0  0.5  0.5    0
#[3,]    0    0  0.5  0.5    0

# Clusters with value within (0,1) are simplified, notice that
# the 0.5s are turned into 1s for only one of the pair. 
#> clust
#     [,1] [,2] [,3] [,4] [,5]
#[1,]    1    1    0    0    1
#[2,]    0    0    1    1    0
#[3,]    0    0    0    0    0

# Use row indices to differentiate clusters by changing non-0s to
# corresponding row indices
#> clust
#      [,1] [,2] [,3] [,4] [,5]
#[1,]    1    1    0    0    1
#[2,]    0    0    2    2    0
#[3,]    0    0    0    0    0

# Cluster grouping of each node is based on colSums of clust directly above
#> clust.grp
#[1] 1 1 2 2 1 # 1-2-5 denoted by 1; 3-4 denoted by 2

