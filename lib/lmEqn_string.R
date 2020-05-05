################################################################################
# GET EQUATION AND R-SQUARED AS STRING
# Modified version of SOURCE: https://groups.google.com/forum/#!topic/ggplot2/1TgH-kG5XMA
################################################################################
lmEqn_string <- function(x="kmer", y="align", data=df){
  m <- eval(parse(text=paste0("lm(", y, "~", x, ", data=data)"
  )
  ))
  p <- cor(x=data[,x], y=data[,y])
  
  a=format(unname(coef(m)[1]), digits=5)
  b=format(unname(coef(m)[2]), digits=5)
  r2=format(summary(m)$r.squared, digits=5)
  p=format(p, digits=5)
  
  #eq <- list(
  #  bquote(bold(bolditalic(y) == .(a) + .(b)~bolditalic(x))),
  #  bquote(bold(bolditalic(r)^2~"="~.(r2)~","~bolditalic(r)~"="~.(p)))
  #)
  eq <- list(bquote(
      atop(italic(y) == .(a) + .(b)~italic(x),
           italic(R)^2~"="~.(r2)~","~italic(R)~"="~.(p)
      )
  ))
  return(eq)
  
  # ORIGINAL
  #eq <- substitute(italic(y) == a + b~italic(x)*","~~italic(r)^2~"="~r2, 
  #                 list(a=format(unname(coef(m)[1]), digits=4),
  #                      b=format(unname(coef(m)[2]), digits=4),
  #                      r2=format(summary(m)$r.squared, digits=3))
  #)
  #as.character(as.expression(eq))
}

