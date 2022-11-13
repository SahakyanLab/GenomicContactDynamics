################################################################################
# Convenient liftOver wrapper as independently implemented (claimed to be with a
# better performance) in the bioconductor package "rtrackplayer"
################################################################################
# REQUIRES: bioconductor package "rtracklayer" and
#           TrantoR function liftOverLoadChainHumanTo from
#           Gen_liftOverLoadChain.R, in case getchain = TRUE.
# BiocManager::install("rtracklayer")

# NOTE
## (1) For one-to-one conversion of coordinates, do liftovers of start and end 
## coordinates separately.
## (2) If ranges are supplied instead, each range can be converted to >1 ranges
## due to insertion in target version or can be shortened due to deletion. 
## Decide based on objective. 
## (3) Liftover reports if a region is inverted or not by inverting the strand,
## but by default, the start and end are already flipped for "-" strands such that
## start < end. If there is no strand information, indicating "+" or "*" for all
## strands will give you same number of regions.
## (4) Output is one-based becacuse start and end can be the same coordinate 
## when referring to a base (width=1). 
## UPDATE: Output will match the coordinate system of the input
## but width reported will always be (end - start + 1) which suggests that 
## inputting one-based coordinates is better and less confusing.
## (5) Note that seqnames and strand are returned as factors.
################################################################################
liftOveR <- function(conversion = "hg19ToHg38", # "hg19ToHg38", "hg19ToDanRer10"
                     space = c("chr2", "chr2", "chr1", "chr3"), # same as chrom
                     start = 708000:708003,
                     end   = 708001:708004,
                     strand = c("+","-","-","+"),
                     getchain = TRUE,
                     rmchain  = TRUE,
                     returnGRangesList = TRUE){

  if(getchain==TRUE){ liftOverLoadChain(chainname=conversion) }

  chain <- import.chain(paste0(conversion,".over.chain"))

  if(rmchain==TRUE){ file.remove(paste0(conversion,".over.chain")) }

  query.granges <- GRanges(seqnames = Rle( space ),
                           ranges = IRanges(start = start, end = end),
                           strand = Rle( strand ) )

  matching.granges <- liftOver(query.granges, chain) # GRangesList object

  if(returnGRangesList==TRUE){
    return(matching.granges)
  } else {
    return(as.data.frame(matching.granges)[,c(1,3:7)]) # returns a simpler dataframe
    # The query entries that did not match are simply absent, hence the column
    # <group> should be used for the correct matching.
    # Example:
    #    group seqnames  start    end width strand
    #  1     1     chr2 708000 708000     1      +
    #  2     2     chr2 708001 708001     1      -
    #  3     3     chr1 772622 772622     1      -
    #  4     4     chr3 666319 666319     1      +
  }

}

################################################################################
# Example of available liftover files for hg19 from:
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/liftOver/
# 
# hg19ToAilMel1.over.chain.gz  04-Feb-2010 16:31   78M
# hg19ToAllMis1.over.chain.gz  28-Jun-2013 14:03   17M
# hg19ToAnoCar1.over.chain.gz  30-May-2009 22:56  8.0M
# hg19ToAnoCar2.over.chain.gz  19-Apr-2011 13:16  7.8M
# hg19ToBalAcu1.over.chain.gz  19-Aug-2016 03:03   80M
# hg19ToBosTau4.over.chain.gz  04-Jun-2009 12:27   81M
# hg19ToBosTau6.over.chain.gz  16-May-2011 20:08   82M
# hg19ToBosTau7.over.chain.gz  23-Jan-2012 14:42   81M
# hg19ToCalJac1.over.chain.gz  14-May-2009 12:57   52M
# hg19ToCalJac3.over.chain.gz  11-Feb-2010 13:46   50M
# hg19ToCanFam2.over.chain.gz  14-May-2009 01:08   88M
# hg19ToCanFam3.over.chain.gz  04-Jul-2012 02:20   87M
# hg19ToCavPor3.over.chain.gz  04-Jun-2009 15:20   83M
# hg19ToCerSim1.over.chain.gz  19-Oct-2012 14:06   82M
# hg19ToChoHof1.over.chain.gz  05-Jun-2009 14:20   67M
# hg19ToCriGri1.over.chain.gz  01-Sep-2013 02:51   71M
# hg19ToDanRer5.over.chain.gz  26-May-2009 19:49  6.4M
# hg19ToDanRer6.over.chain.gz  10-Jul-2009 11:40  7.7M
# hg19ToDanRer7.over.chain.gz  18-Dec-2010 12:42  6.9M
# hg19ToDanRer10.over.chain.gz 25-Sep-2015 08:30  7.2M
# hg19ToDasNov2.over.chain.gz  27-May-2009 14:12   67M
# hg19ToDasNov3.over.chain.gz  09-Jul-2013 13:45   83M
# hg19ToDipOrd1.over.chain.gz  29-May-2009 19:38   58M
# hg19ToEchTel1.over.chain.gz  29-Jun-2012 16:16   54M
# hg19ToEchTel2.over.chain.gz  17-Jun-2013 22:34   61M
# hg19ToEquCab2.over.chain.gz  04-Jun-2009 13:02   82M
# hg19ToEriEur1.over.chain.gz  30-May-2009 07:11   46M
# hg19ToEriEur2.over.chain.gz  09-Jul-2013 08:24   53M
# hg19ToFelCat3.over.chain.gz  04-Jun-2009 15:32   67M
# hg19ToFelCat4.over.chain.gz  08-Jun-2010 11:22   75M
# hg19ToFelCat5.over.chain.gz  27-Feb-2014 13:53   86M
# hg19ToFr2.over.chain.gz      20-May-2009 16:19  3.9M
# hg19ToGalGal3.over.chain.gz  14-May-2009 09:12  7.4M
# hg19ToGalGal4.over.chain.gz  29-Jun-2013 23:12  7.5M
# hg19ToGasAcu1.over.chain.gz  14-May-2009 13:56  4.5M
# hg19ToGeoFor1.over.chain.gz  29-Jul-2012 14:41  7.2M
# hg19ToGorGor1.over.chain.gz  22-Mar-2009 02:59   27M
# hg19ToGorGor3.over.chain.gz  17-Oct-2011 12:44   15M
# hg19ToHg17.over.chain.gz     09-Nov-2012 15:45  341K
# hg19ToHg18.over.chain.gz     04-Jun-2009 15:37  221K
# hg19ToHg38.over.chain.gz     31-Dec-2013 23:08  222K
# hg19ToLoxAfr3.over.chain.gz  22-Jul-2009 12:24   79M
# hg19ToMacEug1.over.chain.gz  31-May-2009 04:35   18M
# hg19ToMacFas5.over.chain.gz  30-Jun-2013 18:50   35M
# hg19ToMelGal1.over.chain.gz  28-Mar-2011 13:01  5.4M
# hg19ToMicMur1.over.chain.gz  22-May-2009 19:16   71M
# hg19ToMm9.over.chain.gz      13-May-2009 20:56   75M
# hg19ToMm10.over.chain.gz     07-Mar-2012 18:40   75M
# hg19ToMonDom5.over.chain.gz  29-May-2009 06:34   32M
# hg19ToMyoLuc1.over.chain.gz  04-Jun-2009 18:38   63M
# hg19ToMyoLuc2.over.chain.gz  30-Aug-2013 01:44   69M
# hg19ToNomLeu1.over.chain.gz  05-Nov-2011 01:05   26M
# hg19ToNomLeu3.over.chain.gz  22-Mar-2013 16:40   26M
# hg19ToOchPri2.over.chain.gz  30-May-2009 00:54   59M
# hg19ToOchPri3.over.chain.gz  18-Jun-2013 03:01   67M
# hg19ToOrnAna1.over.chain.gz  27-May-2009 01:17   18M
# hg19ToOryCun1.over.chain.gz  29-May-2009 06:43   67M
# hg19ToOryCun2.over.chain.gz  12-Aug-2009 21:48   78M
# hg19ToOryLat2.over.chain.gz  22-May-2009 18:25  4.4M
# hg19ToOtoGar1.over.chain.gz  15-May-2009 02:08   77M
# hg19ToOviAri1.over.chain.gz  16-Apr-2010 16:21   57M
# hg19ToOviAri3.over.chain.gz  26-Jun-2013 14:32   81M
# hg19ToPanTro2.over.chain.gz  19-Mar-2009 23:34   13M
# hg19ToPanTro3.over.chain.gz  22-Feb-2011 16:33   13M
# hg19ToPanTro4.over.chain.gz  25-Jan-2013 14:29   13M
# hg19ToPapAnu2.over.chain.gz  19-Aug-2013 23:49   35M
# hg19ToPapHam1.over.chain.gz  20-May-2009 17:56   40M
# hg19ToPetMar1.over.chain.gz  14-May-2009 11:24  3.0M
# hg19ToPetMar2.over.chain.gz  17-Oct-2012 12:59  2.8M
# hg19ToPonAbe2.over.chain.gz  14-May-2009 03:22   26M
# hg19ToProCap1.over.chain.gz  26-May-2009 12:27   64M
# hg19ToPteVam1.over.chain.gz  04-Jun-2009 09:50   75M
# hg19ToRheMac2.over.chain.gz  14-May-2009 08:53   34M
# hg19ToRheMac3.over.chain.gz  15-Mar-2012 19:48   36M
# hg19ToRn4.over.chain.gz      13-May-2009 20:44   70M
# hg19ToRn5.over.chain.gz      27-Jun-2012 17:52   67M
# hg19ToRn6.over.chain.gz      08-Jun-2015 13:42   67M
# hg19ToSaiBol1.over.chain.gz  30-Jun-2013 17:42   48M
# hg19ToSarHar1.over.chain.gz  10-Jul-2013 13:00   16M
# hg19ToSorAra1.over.chain.gz  29-May-2009 08:40   46M
# hg19ToSorAra2.over.chain.gz  18-Jun-2013 16:24   58M
# hg19ToSpeTri1.over.chain.gz  04-Jun-2009 20:15   69M
# hg19ToSpeTri2.over.chain.gz  28-Aug-2013 23:27   87M
# hg19ToSusScr2.over.chain.gz  27-Mar-2010 02:13   71M
# hg19ToSusScr3.over.chain.gz  29-Aug-2013 00:59   79M
# hg19ToTaeGut1.over.chain.gz  26-May-2009 13:04  7.2M
# hg19ToTaeGut2.over.chain.gz  17-Jun-2013 23:33   13M
# hg19ToTarSyr1.over.chain.gz  15-May-2009 16:52   78M
# hg19ToTetNig1.over.chain.gz  14-May-2009 13:41  4.9M
# hg19ToTetNig2.over.chain.gz  10-Aug-2009 14:17  4.0M
# hg19ToTupBel1.over.chain.gz  01-Jun-2009 11:22   69M
# hg19ToTurTru1.over.chain.gz  04-Jun-2009 13:51   77M
# hg19ToVicPac1.over.chain.gz  04-Jun-2009 22:02   67M
# hg19ToXenTro2.over.chain.gz  27-May-2009 11:15  7.5M
# hg19ToXenTro3.over.chain.gz  20-Sep-2011 17:11  7.1M


