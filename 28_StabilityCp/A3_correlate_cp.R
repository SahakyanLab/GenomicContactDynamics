################################################################################
# Correlate recalculated scaled Cp with original scaled Cp
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
wk.dir = file.path(home.dir, "git/GenomicContactDynamics/28_StabilityCp")
recalc_cp.dir = file.path(wk.dir, "out_leave_one_out_cp")
leave_out_list_path = file.path(wk.dir, "out_leave_out_list/leave_out_list.rds")
out.dir = file.path(wk.dir, "out_correlate_cp")
### OTHER SETTINGS #############################################################
gcb = "min2Mb"
chrs = paste0("chr", c("X", 1:22))
plot_only = TRUE
regenerate_corr = TRUE

leave_out_lst <- readRDS(leave_out_list_path)
drop_ids <- names(leave_out_lst)
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggsci)

source(file.path(lib, "doCorTest.R"))
source(file.path(lib, "GG_bgr.R"))
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
if (!plot_only) {
  
  orig_cp_chrcombined <- unlist(
    lapply(chrs, function(chr) {
      readRDS(file.path(recalc_cp.dir, paste(gcb, chr, "original_cp.rds", sep = "_")))
    })
  )
  
  cor_coef_mx <- matrix(NA, nrow = length(drop_ids), ncol = 2,
                        dimnames = list(drop_ids, c("pear", "spea")))
  cor_pval_mx <- matrix(NA, nrow = length(drop_ids), ncol = 2, 
                        dimnames = list(drop_ids, c("pear", "spea")))
  
  for (drop_id in drop_ids) {
    
    recalc_cp_chrcombined <- unlist(
      lapply(chrs, function(chr) {
        readRDS(file.path(recalc_cp.dir, paste(drop_id, gcb, chr, "recalculated_cp.rds", sep = "_")))
      })
    )
    
    # Use package with p-value
    
    if (regenerate_corr) {
      cor_out <- doCorTest(orig_cp_chrcombined, recalc_cp_chrcombined, alt = "two.sided",
                           out.dir = out.dir, out.name = paste0(drop_id, "_", gcb), 
                           method.names = c("pearson", "spearman"))
    } else {
      load(file.path(out.dir, paste0(drop_id, "_", gcb, "_cortest.RData")))
      cor_out <- TEST; rm(TEST)
    }
    
    for (cor_meth in c("pear", "spea")) {
      
      cor_coef_mx[drop_id, cor_meth] <- unname(cor_out[[cor_meth]]$estimate)
      cor_pval_mx[drop_id, cor_meth] <- unname(cor_out[[cor_meth]]$p.value)
      
    }
    
    message(gcb, " ", drop_id, ": done...")
    
  }
  
  # Adjust p-values
  
  cor_padj_mx <- apply(cor_pval_mx, MARGIN = 2, p.adjust, method = "BH")
  
  # Plot correlation values
  
  ## Prepare data for plotting
  cor_coef_df <- cor_coef_mx %>% 
    as.data.frame() %>% 
    rownames_to_column("drop_id") %>% 
    pivot_longer(-drop_id, names_to = "cor_method", values_to = "coef") %>% 
    unite("drop_id_cor_method", c("drop_id", "cor_method"), sep = "_", remove = FALSE)
  
  cor_padj_df <- cor_padj_mx %>% 
    as.data.frame() %>% 
    rownames_to_column("drop_id") %>% 
    pivot_longer(-drop_id, names_to = "cor_method", values_to = "padj") %>% 
    unite("drop_id_cor_method", c("drop_id", "cor_method"), sep = "_", remove = FALSE)
  
  if (identical(cor_coef_df$drop_id_cor_method, cor_padj_df$drop_id_cor_method)) {
    cor_df <- cbind(cor_coef_df, padj = cor_padj_df$padj)
  }
  saveRDS(cor_df, file.path(out.dir, "correlation.rds"))
  
} else {
  
  cor_df <- readRDS(file.path(out.dir, "correlation.rds"))

}

# Plot

num_dropped_celltypes <- lengths(leave_out_lst)
cor_df$drop_n <- num_dropped_celltypes[cor_df$drop_id]
cor_df$drop_n <- factor(as.character(cor_df$drop_n), levels = sort(unique(cor_df$drop_n)))
cor_df$group <- unlist(lapply(strsplit(cor_df$drop_id, "-"), function(x) x[1]))

p <- ggplot(cor_df, aes(x = coef)) +
  geom_density(aes(fill = cor_method, colour = cor_method), alpha = 0.6) +
  scale_colour_jama() +
  scale_fill_jama() +
  bgr2
ggsave(file.path(out.dir, "density_correlation.pdf"), plot = p, height = 10, width = 10)

p <- ggplot(cor_df, aes(x = reorder(group, -coef), y = coef)) +
  geom_violin(scale = "width", trim = TRUE, width = 0.5, colour = "gray80", fill = "gray80") +
  geom_jitter(aes(fill = drop_n, colour = padj < 0.0001), width = 0.1, alpha = 0.7, size = 4, shape = 21) +
  #scale_fill_jama() +
  scale_fill_brewer(palette = "Spectral", direction = -1) +
  scale_colour_manual(values = c("black", "white")) +
  scale_y_continuous(limits = c(round(min(cor_df$coef), 1), 1)) +
  facet_wrap(~ cor_method, ncol = 2) +
  bgr2 +
  theme(axis.text.x = element_text(angle = 45, size = 10, hjust = 1))
ggsave(file.path(out.dir, "violin_correlation.pdf"), plot = p, height = 10, width = 15)

# rm(list=ls()); gc()


# leave_out_lst_1 <- lapply(leave_out_lst[-(1:21)], function(x) {
#   union(x, leave_out_lst$drop_abbrs_outliers_1nmads)
# })

