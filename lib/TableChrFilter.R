################################################################################
#Select chromosomes
################################################################################
TableChrFilter <- function(table = "Robj", chrCol = "Chrcolumnname"){
  chrom.name <- as.character(table[,chrCol])
  table <- table[which( chrom.name == "chr1"  |
                        chrom.name == "chr2"  |
                        chrom.name == "chr3"  |
                        chrom.name == "chr4"  |
                        chrom.name == "chr5"  |
                        chrom.name == "chr6"  |
                        chrom.name == "chr7"  |
                        chrom.name == "chr8"  |
                        chrom.name == "chr9"  |
                        chrom.name == "chr10" |
                        chrom.name == "chr11" |
                        chrom.name == "chr12" |
                        chrom.name == "chr13" |
                        chrom.name == "chr14" |
                        chrom.name == "chr15" |
                        chrom.name == "chr16" |
                        chrom.name == "chr17" |
                        chrom.name == "chr18" |
                        chrom.name == "chr19" |
                        chrom.name == "chr20" |
                        chrom.name == "chr21" |
                        chrom.name == "chr22" |
                        chrom.name == "chrX"  |
                        chrom.name == "chrY"    ), ]
  return(table)
}
################################################################################