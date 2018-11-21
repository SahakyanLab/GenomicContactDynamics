

# > mx <- matrix(1:9, ncol=3, nrow=3)
# > mx
#      [,1] [,2] [,3]
# [1,]    1    4    7
# [2,]    2    5    8
# [3,]    3    6    9
# > library(reshape2)
# > melt(mx)
#   Var1 Var2 value
# 1    1    1     1
# 2    2    1     2
# 3    3    1     3
# 4    1    2     4
# 5    2    2     5
# 6    3    2     6
# 7    1    3     7
# 8    2    3     8
# 9    3    3     9
#
# > acast(melt(mx), Var1~Var2+value)
# 1_1 1_2 1_3 2_4 2_5 2_6 3_7 3_8 3_9
# 1   1  NA  NA   4  NA  NA   7  NA  NA
# 2  NA   2  NA  NA   5  NA  NA   8  NA
# 3  NA  NA   3  NA  NA   6  NA  NA   9
# > acast(melt(mx), Var1~Var2)
#   1 2 3
# 1 1 4 7
# 2 2 5 8
# 3 3 6 9
# > mx2 <- acast(melt(mx), Var1~Var2)
# > mx2
#   1 2 3
# 1 1 4 7
# 2 2 5 8
# 3 3 6 9
# > is.matrix(mx2)
# [1] TRUE
# > mx2[1,3]
# [1] 7
# > melt(mx2)
# Var1 Var2 value
# 1    1    1     1
# 2    2    1     2
# 3    3    1     3
# 4    1    2     4
# 5    2    2     5
# 6    3    2     6
# 7    1    3     7
# 8    2    3     8
# 9    3    3     9
# > acast(melt(mx), Var2~Var1)
#   1 2 3
# 1 1 2 3
# 2 4 5 6
# 3 7 8 9
#
# > melt.mx <- melt(mx)
# > melt.mx
# Var1 Var2 value
# 1    1    1     1
# 2    2    1     2
# 3    3    1     3
# 4    1    2     4
# 5    2    2     5
# 6    3    2     6
# 7    1    3     7
# 8    2    3     8
# 9    3    3     9
#
# > acast(melt(mx), Var1~Var2)
#   1 2 3
# 1 1 4 7
# 2 2 5 8
# 3 3 6 9
# > acast(melt(mx)[sample(1:9, replace=FALSE),], Var1~Var2)
#   1 2 3
# 1 1 4 7
# 2 2 5 8
# 3 3 6 9
# > acast(melt(mx)[sample(1:9, replace=FALSE),], Var1~Var2)
#   1 2 3
# 1 1 4 7
# 2 2 5 8
# 3 3 6 9































#    Var1 Var2 value
# 2     2    1    NA
# 3     3    1    NA
# 4     4    1    NA
# 5     1    2     0
# 6     2    2    NA
# 7     3    2    NA
# 8     4    2    NA
# 9     1    3     0
# 10    2    3     0
# 11    3    3    NA
# 12    4    3    NA
# 13    1    4     0
# 14    2    4     0
# 15    3    4     0
# 16    4    4    NA
#
#
#
#
#
# matrix.filepath <- "AD.nor.chr1.mat" #"./AD.40Kb.raw.chr1.mat"
# #mx <- read.table(matrix.filepath, header=FALSE, as.is=TRUE)
# mx <- data.table::fread(matrix.filepath, sep="\t", header=F)
#
# mx <- as.matrix(mx)
# dimnames(mx)[[2]] <- NULL
#
#
#
#
# library(reshape2)
# library(ggplot2)
# testData=mx[1:100,1:100]
#
# testData=mx[100:105,]
# ggplot(melt(testData), aes(Var2*40000,Var1*40000, fill=value))+
#        xlab("MHz") + ylab("Threshold") +
#        geom_raster() + scale_fill_gradient(low="#FFFFFF", high="#000000")
#
#
#
# library(reshape2)
# library(ggplot2)
# tdm <- melt(testData)
#
# ggplot(tdm, aes(x = Var2, y = Var1, fill = factor(value))) +
#     labs(x = "MHz", y = "Threshold", fill = "Value") +
#     geom_raster() +
#     scale_fill_manual(breaks = levels(factor(tdm$value)),
#                       values = c("white", "black")) +
#     theme(plot.background = element_rect(fill = "grey90"),
#           legend.background = element_rect(fill = "grey90")) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0))
#
# It was necessary to create a different plot background color because
# you wanted to remove the padding, which blended much of the plot
# boundary with the original white background, and had to do something
# similar to the legend background so that you could see the white
# legend box.
#
# The expand = c(0, 0) argument to the two scale functions answers your
# second question, and converting value from numeric to vector answers
# the first.
