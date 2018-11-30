################################################################################
# This small function downloads and gunzips the chain file based on the supplied
# chain name. Data are retrieved from the UCSC repositories. The way system
# gunzip is used now, makes it only Linux/Unix/MacOS compatible. Switch to
# gunzip from the R.utils package for the Windows compatibility.
################################################################################

liftOverLoadChain <- function(chainname="hg19ToDanRer10"){

  download.file(url=paste0("http://hgdownload.cse.ucsc.edu/goldenPath/",
                           strsplit(chainname,"To")[[1]][1],"/liftOver/",
                           chainname,".over.chain.gz"),
                destfile=paste0("./",chainname,".over.chain.gz"))

  system(paste0("gunzip ./",chainname,".over.chain.gz") )

}

################################################################################
