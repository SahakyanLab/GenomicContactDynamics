# Check distribution of satellite DNA between contacting regions

# Setup

lib = "../lib"
library(assertthat)
library(tidyverse)
library(ggrepel)
source(paste0(lib, "/GG_bgr.R"))

# Parameters

sumrep_dir = "z_ignore_git/out_HicRepeatExploration/subfamSatellite_sumrep"
metric_dir_list = list(
  minrep = "z_ignore_git/out_HicRepeatExploration/subfamSatellite_minrep",
  skewrep = "z_ignore_git/out_HicRepeatExploration/subfamSatellite_skewrep"
)

#chr = "chr21"
chrs = paste0("chr", c(1:22, "X"))
gcb = "min2Mb"
sumrep_min = 2 # For sumrep, consider only contacts with 2 sites in total

persistent_cps = 19:21
dynamic_cps = 1:3
n_samples = 10000
seed_val = 216

metric = "minrep"
out_dir = file.path("z_ignore_git/out_satellite", metric)
dir.create(out_dir, recursive = TRUE)

# ----- MAIN -----

# Get number of contacts with Cp = {1:3, 19:21} (across chromosomes) per repeat subfamily names with at least 2 sites
# i.e. are there enough contacts to compare dynamic vs. persistent contact per repeat subfamily?

for (i in 1:length(chrs)) {
  
  load(file.path(sumrep_dir, paste0(chrs[i], "_MinElm_", gcb, ".RData")))
  sumrep_orig_mx <- MINELM.MX
  sumrep_mx <- sumrep_orig_mx[sumrep_orig_mx[, "ntis"] %in% c(1:3, 19:21), ]
  element_sums <- colSums((sumrep_mx[, -1]) >= sumrep_min)

  if (i == 1){
    collected_sums_tibble <- element_sums
  } else {
    collected_sums_tibble <- bind_rows(collected_sums_tibble, element_sums)
  }
  
}

all_element_sums <- colSums(collected_sums_tibble, na.rm = TRUE)
p <- all_element_sums %>% 
  enframe() %>% 
  ggplot(aes(name, value)) + 
  geom_col() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

# Get per element list with persistent and dynamic metric values

element_names <- names(all_element_sums[all_element_sums > 0])#[1:5] # REMOVE
writeLines(element_names, file.path(out_dir, "integer_element_filename_mapping.txt"))

for (e in 1:length(element_names)) {
  element <- element_names[e]
  values_list <- list(list(persistent = NULL, dynamic = NULL))
  names(values_list) <- element
  for (i in 1:length(chrs)) {
    # Choose contacts with at least 2 sites based on sumrep
    load(file.path(sumrep_dir, paste0(chrs[i], "_MinElm_", gcb, ".RData")))
    sumrep_mx <- MINELM.MX
    if (!element %in% colnames(sumrep_mx[, -1])) {
      next
    } else {
      is_persistent <- (sumrep_mx[, element] >= sumrep_min) & (sumrep_mx[, "ntis"] %in% persistent_cps)
      is_dynamic <- (sumrep_mx[, element] >= sumrep_min) & (sumrep_mx[, "ntis"] %in% dynamic_cps)
      #if (any(c(sum(is_persistent), sum(is_dynamic)) == 0)) { next }
      # metric
      load(file.path(metric_dir_list[[metric]], paste0(chrs[i], "_MinElm_", gcb, ".RData")))
      minrep_mx <- MINELM.MX
      values_list[[element]] <- list(
        persistent = c(values_list[[element]]$persistent, minrep_mx[is_persistent, element]),
        dynamic = c(values_list[[element]]$dynamic, minrep_mx[is_dynamic, element])
      )
    }
  } # chrs for loop
  saveRDS(values_list, file.path(out_dir, paste0("values_", e, ".rds")))
} # element_names for loop

# Get data frame (combining elements) with persistent and dynamic (means of samples) metric values

data_list <- list()
n_points_list <- list()
for (e in 1:length(element_names)) {
#for (e in 1:2) {
  element <- element_names[e]
  values_list <- readRDS(file.path(out_dir, paste0("values_", e, ".rds")))
  n_points <- unlist(lapply(values_list[[element]], length))
  n_points_list[[element]] <- data.frame(element, t(data.frame(rev(n_points))))
  message(element, "...; n_points:")
  print(n_points)
  if (all(n_points == 0)) {
    next; message("Skipping ", element)
  } else {
    if (n_points[["persistent"]] >= 5) {
      # Sample
      set.seed(seed_val)
      mean_of_samples <- sapply(1:n_samples, simplify = FALSE, function(i) {
        sample(
          values_list[[element]]$dynamic, size = length(values_list[[element]]$persistent), replace = FALSE
        ) %>% mean()
      }) %>% unlist()
      values_list[[element]][["dynamic.mean_of_samples"]] <- mean_of_samples
    } else {
      values_list[[element]][["dynamic.mean_of_samples"]] <- NULL
    }
    saveRDS(values_list, file.path(out_dir, paste0("values_withsampleddata_", e, ".rds")))
    
    p_df_list <- list()
    if (length(values_list[[element]]$persistent) > 0) {
      p_df_list[[1]] <- data.frame(cp_type = "persistent", value = values_list[[element]]$persistent)
    }
    if (length(values_list[[element]]$dynamic) > 0) {
      p_df_list[[2]] <- data.frame(cp_type = "dynamic", value = values_list[[element]]$dynamic)
    }
    if (!is.null(values_list[[element]][["dynamic.mean_of_samples"]])) {
      p_df_list[[3]] <- data.frame(cp_type = "dynamic.mean_of_samples", value = values_list[[element]]$dynamic.mean_of_samples)
    }
    p_df <- do.call("rbind", p_df_list)
    rm(p_df_list)
  
    #p_df <- rbind(
    #  data.frame(cp_type = "persistent", value = values_list[[element]]$persistent),
    #  data.frame(cp_type = "dynamic.mean_of_samples", value = values_list[[element]]$dynamic.mean_of_samples)
    #)
    p_df <- cbind(element = element, p_df)
    data_list[[element]] <- p_df
  }
}

#
#p1_df <- data_list$TAR1
#p1_df %>% 
#  filter(cp_type == "dynamic.mean_of_samples") %>% 
#  pull(value) %>% 
#  mean()

n_points_df <- do.call("rbind", n_points_list)
data_df <- do.call("rbind", data_list)
# Plot distribution of values in values_list
p_list <- list()
pdf(file.path(out_dir, "distribution.pdf"), height = 5, width = 5)
par(mfrow = c(1, 1))
for (element_lvl in unique(data_df$element)) {
  p_list[[element_lvl]] <- data_df %>% 
    filter(element == element_lvl) %>% 
    ggplot(aes(value)) +
    geom_bar(colour = "black") +
    #facet_grid(element ~ .) +
    facet_grid(cp_type ~ .) +
    theme_classic()
  print(p_list[[element_lvl]])
}
dev.off()

# Create data for heatmap showing mean of metric per element for persistent and dynamic (means of samples)
data_df$cp_type <- factor(
  data_df$cp_type, levels = c("dynamic", "dynamic.mean_of_samples", "persistent")
)
data_agg_df <- data_df %>% 
  group_by(element, cp_type) %>% 
  summarise(
    mean = mean(value),
    median = median(value)
  ) %>% 
  ungroup() %>% 
  arrange(cp_type) %>% 
  pivot_wider(
    names_from = cp_type, values_from = c(median, mean), id_cols = element, names_glue = "{.value}.{cp_type}"
  )
#colnames(data_agg_df) <- c("element", paste0("mean.", colnames(data_agg_df)[-1]))
colnames(n_points_df) <- c("element", paste0("n.", colnames(n_points_df)[-1]))
combined_df <- n_points_df %>% 
  left_join(data_agg_df, join_by(element))

# Significance - one sided observed > null

elements_with_enough_datapoints <- combined_df %>% 
  filter(n.persistent >= 5) %>% 
  pull(element)
pvals <- sapply(elements_with_enough_datapoints, function(element_lvl) {
  element_df <- data_df %>% 
    filter(element == element_lvl)
  # How many null values > observed
  mean_observed <- combined_df[combined_df$element == element_lvl, "mean.persistent"]
  null_values <- element_df %>% 
    filter(cp_type == "dynamic.mean_of_samples") %>% 
    pull(value)
  assert_that(
    length(null_values) == n_samples,
    msg = paste0(element_lvl, ": Error length(null_values) == n_samples")
  )
  if (metric == "minrep") {
    pval <- sum(null_values > mean_observed) / length(null_values)
  } else if (metric == "skewrep") {
    pval <- sum(null_values < mean_observed) / length(null_values)
  } else {
    stop("Unexpected metric argument")
  }
  return(pval)
})
pvals_df <- pvals %>% 
  enframe(name = "element", value = paste0("pvalue_basedon_", n_samples, "permutedsamples"))
padjs <- pvals_df[[paste0("pvalue_basedon_", n_samples, "permutedsamples")]]
pvals_df$padj <- p.adjust(padjs, method = "BH")

combined_df <- combined_df %>% 
  left_join(pvals_df, join_by(element))
write.csv(combined_df, file.path(out_dir, "combined_stats.csv"), row.names = FALSE)

# rm(list = ls()) ; gc()

# Plot

df.plot <- read_csv(file.path(out_dir, "combined_stats.csv"))
df.plot$label <- df.plot$element
withLabel = TRUE
#df.plot <- DF[DF$plot.group %in% c("ijcov", "aminrepnon0.fr"),]
#df.plot$label <- NA
#is.label <- !is.na(df.plot$per) & df.plot$per > 0.1
#df.plot$label[is.label] <- as.character(df.plot$repName)[is.label]
#p <- ggplot(data=df.plot, aes(x=dyn, y=per)) +
p <- ggplot(data=df.plot, aes(x=mean.dynamic, y=mean.persistent)) +
  #geom_hline(yintercept=0.5, linetype="dashed") +
  #geom_vline(xintercept=0.5, linetype="dashed") +
  geom_abline(slope=1, intercept=0, linetype="dashed", colour="black") +
  #geom_hex() +
  #stat_binhex(aes(label=..count..), geom="text", colour="white") +
  geom_point(shape=20, size=10, col=adjustcolor("black", 0.1)) +
  geom_point(shape=20, size=2, col=adjustcolor("black", 0.5)) +
  scale_x_continuous(limits=c(0,1.0)) +
  scale_y_continuous(limits=c(0,1.0)) +
  #scale_fill_gradient(low="#3C5488FF", high="#DC0000FF") + 
  bgr2 #+
  #facet_grid(.~plot.group)

if(withLabel){
  p <- p +
    geom_text_repel(aes(label=label), box.padding=1, col="darkgreen",
                    size=1, segment.size=1, min.segment.length=0) 
} else {
  p <- p + 
    theme(axis.title.y=element_blank(), axis.text.x=element_blank(), 
          axis.text.y=element_blank())
}

out.affix = paste0(metric, "_", withLabel)
out.dir = out_dir
ggsave(filename=paste0(out.dir, "/scatter_dynCp1To3VsPerCp19To21_", out.affix, ".png"),
       plot=p, height=10*300, width=20*300, units="px")
