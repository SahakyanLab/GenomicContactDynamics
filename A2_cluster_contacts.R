# Here we shall load all the cross-tissue/cell contact data and cluster all the
# contacts as per their

contact.path = "/Volumes/Data/Database/GSE87112/combined_contacts/HiCNorm_QQ_primary_cohort"
chr = 1
nCPU = 4





load(paste0(contact.path,"/human_chr",chr,"_allcontacts.RData"))



MELT.MX$upper.tri[,-c(1,2)]

library(cluster)




#MELT.MX$upper.tri[,-c(1,2)][1:10, 1:21]


row.means   <- rowMeans(MELT.MX$upper.tri[,-c(1,2)])
row.sds     <- apply(MELT.MX$upper.tri[,-c(1,2)], 1, sd)
row.medians <- apply(MELT.MX$upper.tri[,-c(1,2)], 1, median)
row.num0s   <- apply(MELT.MX$upper.tri[,-c(1,2)], 1, FUN=function(i){sum(i==0)})



plot(x=row.means[row.num0s<1], y=row.sds[row.num0s<1], cex=0.01, xlim=c(0,10), ylim=c(0,10))
plot(x=row.means[row.num0s<=10], y=row.sds[row.num0s<=10], cex=0.01, xlim=c(0,10), ylim=c(0,10))

plot(x=row.means[row.num0s<=5], y=row.sds[row.num0s<=5], cex=0.1)






MELT.MX$upper.tri[,-c(1,2)]


require(parallel)
res <- mclapply( 1:ncol(origmatrix) , mc.cores = 1 ,
                   function(x){ c( mean( origmatrix[,x] ) ,
                   sd( origmatrix[,x] ) , var( origmatrix[,x] ) ) } )

# So the first element of the resulting list looks like
 res[[1]]
  # [1] 5.500000 3.027650 9.166667

df <- as.data.frame( res )
rownames(df) <- c("mean","sd","var")
colnames(df) <- colnames(origmatrix)
#               a         b         c         d         e
#   mean 5.500000 15.500000 25.500000 35.500000 45.500000
#   sd   3.027650  3.027650  3.027650  3.027650  3.027650
#   var  9.166667  9.166667  9.166667  9.166667  9.166667



a <- mclapply(X=1:nrow(MELT.MX$upper.tri),
         FUN=function(x){
           sd( MELT.MX$upper.tri[x,-c(1,2)] )
         },
         mc.preschedule = FALSE, mc.set.seed = FALSE,
         mc.silent = FALSE, mc.cores = nCPU,
         mc.cleanup = TRUE, mc.allow.recursive = TRUE,
         affinity.list = NULL)






mcmapply(FUN, ...,
         MoreArgs = NULL, SIMPLIFY = TRUE, USE.NAMES = TRUE,
         mc.preschedule = TRUE, mc.set.seed = TRUE,
         mc.silent = FALSE, mc.cores = getOption("mc.cores", 2L),
         mc.cleanup = TRUE, affinity.list = NULL)




fraction <- 0.01

set.seed(845)
sample.ind <- sample(x=1:dim(MELT.MX$upper.tri)[1], replace=FALSE,
                     size=floor(dim(MELT.MX$upper.tri)[1]*fraction))


#explore the tSNE visualisation and Clara clustering possibilities



a <- clara(
  x = MELT.MX$upper.tri[,-c(1,2)],
  k = 10,
  metric = "euclidean",
  samples=500,
  sampsize=100,
  stand=FALSE
)





library(Rtsne)


tsne <- Rtsne(MELT.MX$upper.tri[sample.ind,-c(1,2)],
              dims = 2, perplexity=30, verbose=TRUE,
              max_iter = 5000, check_duplicates=FALSE, pca=TRUE)





library(reshape2)
save(MELT.MX, file=paste0(out.path,"/human_chr",chr,"_allcontacts.RData"))
