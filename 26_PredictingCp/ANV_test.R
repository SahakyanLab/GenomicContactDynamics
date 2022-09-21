# Sequence ATCTAGCT
seq = c("A", "T", "C", "T", "A", "G", "C", "T")
AIF <- rbind(
  u_Ai=c(1,1,1,1,2,2,2,2),
  u_Ci=c(0,0,1,1,1,1,2,2),
  u_Gi=c(0,0,0,0,0,1,1,1),
  u_Ti=c(0,1,1,2,2,2,2,3)
)
dimnames(AIF)[[2]] <- seq

u_base=rowMeans(AIF)
n_base <- table(seq)

 
base1 = "A"
base2 = "C"
id.1 <- paste0("u_", base1, "i")
id.2 <- paste0("u_", base2, "i")
tmp <- sum( (AIF[id.1,] - u_base[id.1]) * (AIF[id.2,] - u_base[id.2]) ) 
COVb1b2 <- unname( tmp / (n_base[base1] * n_base[base2]) )
COVb1b2

COV:
AC = 0.5
AG = 0.75
AT = 0.41667
CG = 1
CT = 0.66667
GT = 0.70833
