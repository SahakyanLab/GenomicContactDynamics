################################################################################
# Per gcb and gap, compare complementarity values between persistent and variable
# contacts
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar = TRUE) 
options(warn = 1) 

whorunsit = "LiezelMac" # "LiezelMac", "LiezelCluster", "LiezelLinuxDesk",
# "AlexMac", "AlexCluster"

if( !is.null(whorunsit[1]) ){
  # This can be expanded as needed ...
  if(whorunsit == "LiezelMac"){
    home.dir = "/Users/ltamon"
    os = "Mac"
  } else if(whorunsit == "LiezelCluster"){
    home.dir = "/project/sahakyanlab/ltamon"
    os = "Linux"
  } else {
    stop("The supplied <whorunsit> option is not created in the script.", quote=F)
  }
}
lib = file.path(home.dir, "git/GenomicContactDynamics/lib/")
wk.dir = file.path(home.dir, "git/GenomicContactDynamics/25_LengthDependence")
data.dir = file.path(home.dir, "data/gcd/Database")
src.dir = file.path(wk.dir, "out_get_complementarity_perCpGap")
out.dir = file.path(wk.dir, "out_compare_complementarity_perCpGap")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chrs = paste0("chr", c("X", 22:1)) #paste0("chr", c("X", 1:22))
gaps = 50:1000
compl.type = "align" # kmer | align | Gfree
plot_only = TRUE
min_ij_group = 100 # Compare at gap with >= 100 persistent and >= 100 variable contacts
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(dplyr)
library(tibble)
library(tidyr)
library(stringr)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
source(file.path(lib, "compareTwoDist.R"))
source(file.path(lib, "GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if (!plot_only) {
  
  mx <- matrix(data = NA, nrow = length(gaps), ncol = 11,
               dimnames = list(paste0("gap_", as.character(gaps)),
                               c("num_nonNAcompl_var.Cp", "num_nonNAcompl_per.Cp",
                                 "mean_var.Cp", "mean_per.Cp",
                                 "median_var.Cp", "median_per.Cp",
                                 "sd_var.Cp", "sd_per.Cp",
                                 paste0("pval_", c("t", "mw", "ks")))
               )
  )
  
  for (gap in gaps) {
    
    gap_char <- as.character(gap)
    compl_per_gap_lst <- vector(mode = "list", length = 2)
    names(compl_per_gap_lst) <- c("var.Cp", "per.Cp")
    
    for (chr in chrs){
      
      obj_path <- file.path(src.dir, 
                            paste0(gcb, "_gap", gap_char, "_", chr, "_", compl.type, "_per_cpgroup.rds"))
      if (file.exists(obj_path)) {
        
        compl_per_chr_per_gap_lst <- readRDS(obj_path)
        compl_per_gap_lst[["var.Cp"]] <- c(compl_per_gap_lst[["var.Cp"]], 
                                           compl_per_chr_per_gap_lst[["var.Cp"]])
        compl_per_gap_lst[["per.Cp"]] <- c(compl_per_gap_lst[["per.Cp"]], 
                                           compl_per_chr_per_gap_lst[["per.Cp"]])
        
      }
      
    }
    
    if (all(is.null(unlist(compl_per_gap_lst)))) {
      next
    }
    
    # Per gap, number of persistent and variable contacts used for comparison
    TEST <- tryCatch(
      {
        compareTwoDist(x = compl_per_gap_lst$var.Cp, y = compl_per_gap_lst$per.Cp)  
      },
      error = function(e) {
        TEST <- list(t = list(p.value = NA), mw = list(p.value = NA), ks = list(p.value = NA))
        return(TEST)
      }
    )
    
    mx[paste0("gap_", gap_char), ] <- unname(
      c(vapply(compl_per_gap_lst, function(x) sum(!is.na(x)), FUN.VALUE = numeric(1)),
        vapply(compl_per_gap_lst, mean, FUN.VALUE = numeric(1), na.rm = TRUE),
        vapply(compl_per_gap_lst, median, FUN.VALUE = numeric(1), na.rm = TRUE),
        vapply(compl_per_gap_lst, sd, FUN.VALUE = numeric(1), na.rm = TRUE),
        vapply(c("t", "mw", "ks"), FUN = function(x) TEST[[x]][["p.value"]], FUN.VALUE = numeric(1))
      )
    )
    
    rm(compl_per_gap_lst, TEST)
    message(gcb, " gap=", gap_char, " done!")
    
  }
  
  mx <- mx[!is.na(mx[,1]) & !is.na(mx[,2]),]
  saveRDS(mx, file.path(out.dir, paste0(gcb, "_", compl.type, "_compare_complementarity.rds")))
  
} else {
  mx <- readRDS(file.path(out.dir, paste0(gcb, "_", compl.type, "_compare_complementarity.rds")))
}

compare_mx <- mx[mx[,"num_nonNAcompl_var.Cp"] >= min_ij_group & mx[,"num_nonNAcompl_per.Cp"] >= min_ij_group,]

# Prepare plotting data

# Adjust p-values
for (pval_colname in grep("pval_", colnames(compare_mx), value = TRUE)) {
  compare_mx[ ,pval_colname] <- p.adjust(compare_mx[ ,pval_colname], method = "BH")
}

plot_df <- compare_mx %>% 
  as.data.frame() %>% 
  rownames_to_column("gap_id") %>% 
  pivot_longer(!c(gap_id, num_nonNAcompl_var.Cp, num_nonNAcompl_per.Cp, pval_t, pval_mw, pval_ks)) %>% 
  separate(name, c("metric", "cp_group"), sep = "_", remove = FALSE) %>% 
  #separate(gap_id, c("gap_id", "gap"), sep = "_") %>% 
  mutate(gap_val = as.integer(str_replace(gap_id, "gap_", ""))) %>% 
  filter(!metric %in% c("sd"))
plot_df$gap_group <- cut(plot_df$gap_val, breaks = seq(50, 150, by = 10), 
                         ordered_result = TRUE, include.lowest = TRUE)
plot_df$cp_group <- factor(plot_df$cp_group, levels = c("var.Cp", "per.Cp"))

# Plot

set.seed(290) # For geom_jitter
p <- ggplot(plot_df, aes(x = cp_group, y = value)) +
  geom_boxplot(alpha = 0.5, aes(colour = cp_group)) +  
  geom_jitter(aes(fill = gap_val), width = 0.2, size = 1.5, alpha = 0.7, shape = 21) +
  stat_compare_means(comparisons = list(c("var.Cp", "per.Cp")), ref.group = "var.Cp",
                     method = "wilcox.test", paired = TRUE, label = "p.format", size = 1) +
  scale_colour_manual(values = colorRampPalette(brewer.pal(n = 11, name = "Spectral"))(21)[c(21, 1)]) +
  scale_fill_distiller(palette = "BrBG", direction = 1) +
  labs(title = paste0(gcb, "_", compl.type, ", paired wilcox p-value"), x = NULL) +
  facet_wrap(~ metric, ncol = 2) +
  bgr1
ggsave(file.path(out.dir, paste0(gcb, "_", compl.type, "_boxplot_correlation.pdf")), plot = p, height = 5, width = 15)

# rm(list=ls()); gc()
