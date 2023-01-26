# Load COR R object


# Spearman rho

cor.mx <- lapply(COR, FUN=function(cor.obj){
  c(cor.obj$spea$estimate, cor.obj$spea$p.value)
})
cor.mx <- do.call(rbind, cor.mx)
range(cor.mx[,"rho"])
