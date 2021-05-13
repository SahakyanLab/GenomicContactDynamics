################################################################################
# Plot per location, per filter, combining all signatures
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
start.time <- Sys.time()

# Expands warnings
options(warn=1)

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    lib = "/Users/ltamon/DPhil/lib"
    data.dir = "/Users/ltamon/Database"
    wk.dir = "/Users/ltamon/DPhil/GCD_polished/19_MutationRatesVsPersist"
    os = "Mac"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
plot.dir = paste0(wk.dir, "/out_plotdataVsCp_PCAWG")
out.dir = paste0(wk.dir, "/out_signatureCombinedPlot")

comb1.file = paste0(wk.dir, "/out_makecombfile/combfile_nosampfilter_1000rawInf_-1_0")
comb2.file = paste0(wk.dir, "/out_makecombfile_sigPlot/combfile_min2Mb_donor_centric_PCAWG_Hg19_Nmsitenorm_mean")
numbers.file = paste0(wk.dir, "/out_checkNumbers/donor_centric_PCAWG_Hg19_comprehensivefinal.csv")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
data.id = "donor_centric_PCAWG" # "CosmicNCV", "donor_centric_PCAWG"
src.id = "Hg19" # "Hg19" | "hg38ToHg19"
calc = "Nmsitenorm" # "Tmut" | "Nmsite" | "TmutDIVNmsite" | "Nmsitenorm" | "numWTSEQ"
aggregatefunx = "mean" # "mean" | "median"
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(ggplot2)
library(ggpubr)
library(yarrr)
yarrr::piratepal("basel")
source(paste0(lib, "/GG_bgr.R"))
### FUNCTION ###################################################################
makePlot <- function(df, calc, yint=NULL, plot.id){
  
  if( !is.numeric(df$ind) & !is.character(df$ind) ){
    stop("Plot making: ind not numeric or character as expected.")
  }
  
  if( !is.numeric(df$pval) ){
    stop("Plot making: pval not numeric as expected.")
  }
  
  lev.v <- as.character(sort( as.numeric(unique(df$ind)) ))
  df$ind <- factor(x=as.character(df$ind), levels=lev.v)
  
  df$pval <- c(`0`="ns", `1`="sig")[ as.character(as.numeric(df$pval<0.05)) ]
  df$pval <- factor(x=as.character(df$pval), levels=c("ns", "sig"))
  shape.v <- c(8, 1)
  names(shape.v) <- levels(df$pval)
  
  col.v <- colorRampPalette(yarrr::piratepal("basel"))(length(unique(df$SIG.id)))
  #df$SIG.id <- factor(df$SIG.id, levels=unique(c("nosampfilter", sort(df$SIG.id))))
  # Assign colour black to nosampfilter
  col.v[1] <- "black"
  
  p <- ggplot(data=df, aes(x=ind, y=values, group=SIG.id)) +
    geom_line(aes(colour=SIG.id), alpha=0.3) +
    geom_point(aes(colour=SIG.id, shape=pval), size=2) +
    scale_colour_manual(values=col.v) +
    scale_shape_manual(values=shape.v) + 
    labs(x="Cp", y=calc, title=plot.id) +
    bgr2 + 
    theme(plot.title=element_text(size=15), legend.text=element_text(size=2),
          legend.title=element_text(size=10))
    
  if( !is.null(yint) ){
    p <- p + geom_hline(yintercept=yint, linetype="dashed", colour="gray70", size=0.5) 
  }
  
  return(p)
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
comb.v <- readLines(con=comb1.file)
comb.len <- length(comb.v)

DF <- list()
for(cb in 1:comb.len){
  
  cb.v <- strsplit(x=comb.v[cb], split=";", fixed=T)[[1]]
  names(cb.v) <- c("SIG", "mut.id", "sigEpLim", "loc")
  cb.v["mut.id"] <- gsub(x=cb.v["mut.id"], replacement="To", pattern=">", fixed=T)
  
  pdata.id <- paste0(gcb, "_", aggregatefunx, "_", calc, "_", data.id, "_", 
                     src.id, "_", cb.v["mut.id"], "_", cb.v["SIG"], "_", 
                     cb.v["loc"], "_sigEperclimits_", cb.v["sigEpLim"])
  fle <- paste0(plot.dir, "/", pdata.id, "_scatplot.RData")
  if( file.exists(fle) ){
    
    load(file=fle) # P.DF
    load(file=paste0(plot.dir, "/", pdata.id, "_FC_scatplot.RData")) # PFC.DF
    
  } else {
    
    print(paste0(pdata.id, " skipped!"), quote=F)
    next
    
  }
  
  DF[[pdata.id]] <- cbind.data.frame(PFC.DF, values.abs=P.DF$values,
                                     stringsAsFactors=F)
  
  rm(P.DF, PFC.DF, pdata.id, cb.v, cb)
  
} # comb.len for loop end
rm(comb.v, comb.len); gc()

DF <- do.call("rbind.data.frame", c(DF, stringsAsFactors=F))
rownames(DF) <- NULL

# Separate nosampfilter rows
NS.DF <- DF[DF$sigEpLim.id=="nosampfilter",]
NS.DF$SIG.id <- "nosampfilter"
NS.DF <- NS.DF[!duplicated(NS.DF),]

# Separate nosampfilter data (as NS.DF) from main DF
DF <- DF[DF$sigEpLim.id!="nosampfilter",]

# Get important numbers per MUTBIN.DF ready as legend labels
numbers.df <- read.csv(file=numbers.file, stringsAsFactors=F, header=T)
numbers.df[numbers.df$sigEpLim.id=="nosampfilter", "SIG.id"] <- "nosampfilter"
numbers.df <- numbers.df[!duplicated(numbers.df),]

numbers <- paste0(numbers.df$Nmut, "_", numbers.df$Nsamp)
nsamp <- numbers.df$Nsamp
names(numbers) <- names(nsamp) <- paste0(numbers.df$mut.id, numbers.df$SIG.id, numbers.df$loc.id, numbers.df$sigEpLim.id)

#-------------------Plot
loc.id.v <- unique(DF$loc.id)
loc.id.v.len <- length(loc.id.v)
sigEpLim.id.v <- unique(DF$sigEpLim.id)
sigEpLim.id.v[sigEpLim.id.v!="nosampfilter"]
sigEpLim.id.v.len <- length(sigEpLim.id.v)
mut.id.v <- unique(DF$mut.id)
mut.id.v.len <- length(mut.id.v)

#s.v <- rep(x=sigEpLim.id.v, times=mut.id.v.len)
s.v <- rep(x=loc.id.v, each=mut.id.v.len)
s.v.len <- length(unique(s.v))
#m.v <- rep(x=mut.id.v, each=sigEpLim.id.v.len)
m.v <- rep(x=mut.id.v, times=loc.id.v.len)
len <- unique(c(length(s.v), length(m.v)))

comb.v <- readLines(con=comb2.file)
comb.len <- length(comb.v)

for(cb in 1:comb.len){
  
  cb.v <- strsplit(x=comb.v[cb], split=";", fixed=T)[[1]]
  if(length(cb.v)!=3){ stop("Invalid comb.") }
  wise <- cb.v[1]
  sigEpLim.id <- cb.v[2]
  Nsamp.min <- cb.v[3]

  # Filter: bin-wise or contact-wise?
  df1 <- DF[DF$wise==wise & DF$sigEpLim.id==sigEpLim.id,]
  
  out.id <- paste0(gcb, "_", aggregatefunx, "_", calc, "_", data.id, "_", src.id,
                   "_sigEpLim", sigEpLim.id, "_Nsampmin", Nsamp.min, "_", wise)
  
  #for(loc.id in loc.id.v){
  
  p.lst <- list()
  for(i in 1:len){
    
    mut.id <- m.v[i]
    loc.id <- s.v[i]
    
    df <- df1[df1$loc.id==loc.id & df1$mut.id==mut.id,]
    df <- rbind.data.frame(df, NS.DF[NS.DF$loc.id==loc.id & NS.DF$mut.id==mut.id & NS.DF$wise==wise,])
    df$wise <- NULL
    
    perclim.df <- merge(x=aggregate(x=df$percbinmut, by=list(as.factor(df$SIG.id)), FUN=min, na.rm=F),
                        y=aggregate(x=df$percbinmut, by=list(as.factor(df$SIG.id)), FUN=max, na.rm=F),
                        by="Group.1")
    perclim.df$id <- paste0(mut.id, as.character(perclim.df$Group.1), loc.id, sigEpLim.id)
    perclim.df$numbers <- numbers[perclim.df$id]
    perclim.df$nsamp <- nsamp[perclim.df$id]
    row.names(perclim.df) <- perclim.df$Group.1
    ns.id <- paste0(mut.id, "nosampfilter", loc.id, "nosampfilter")
    perclim.df["nosampfilter","numbers"] <- unname(numbers[ns.id])
    perclim.df["nosampfilter","nsamp"] <- unname(nsamp[ns.id])
    perclim.df$Group.1 <- paste0(perclim.df$Group.1, "-", perclim.df$numbers, "-", 
                                 round(perclim.df$x.x, digits=1), "-", round(perclim.df$x.y, digits=1))
    if( any(is.na(c(perclim.df$numbers, perclim.df$nsamp))) ){
      stop( paste0(mut.id, " ", loc.id, " ", sigEpLim.id, ": Checkpoint 1.") )
    }
    df$nsamp <- perclim.df[df$SIG.id,"nsamp"]
    df$SIG.id <- perclim.df[df$SIG.id,"Group.1"]
    
    # Rows with Nsamp < Nsamp.min should be removed
    df <- df[df$nsamp>=as.numeric(Nsamp.min),]
    
    p.lst[[as.character(i)]] <- makePlot(df=df, calc=calc, yint=0,
                                         plot.id=paste0(out.id, "\n", loc.id, "_", sigEpLim.id, "_", mut.id))
    
    rm(mut.id, loc.id, df, perclim.df, ns.id)
    
    print(i)
    
  }
  
  p.arr <- ggarrange(plotlist=p.lst, nrow=s.v.len, ncol=mut.id.v.len, legend=NULL)
  ggexport(p.arr, height=10*s.v.len, width=10*mut.id.v.len, 
           filename=paste0(out.dir, "/", out.id, ".pdf"))
  
  rm(wise, sigEpLim.id, Nsamp.min, cb.v, p.lst, df1, p.arr, out.id)
  
  #}
  
} # comb.len for loop end

# rm(list=ls()); gc()
