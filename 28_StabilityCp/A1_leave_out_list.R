################################################################################
# R list containing cell type/s to leave out when recalculating Cp
################################################################################
# FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS * FLAGS
### DIRECTORY STRUCTURE ########################################################
options(warnPartialMatchDollar=TRUE) 
options(warn=1)

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
wk.dir = file.path(home.dir, "git/GenomicContactDynamics/28_StabilityCp")
dataqc_path = file.path(wk.dir, "data_quality.csv")
### OTHER SETTINGS #############################################################
cell_info_path = file.path(wk.dir, "celltypes.csv")
celltype_ids = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
                 "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
library(scuttle)

select_represent_group <- function(cell_info_df) {
  
  unique_group_names <- names(which(table(cell_info_df$group) == 1))
  unique_group_abbrs <- cell_info_df$abbr[cell_info_df$group %in% unique_group_names]
  
  # Groups with more than one cell type
  dupl_group_names <- setdiff(cell_info_df$group, unique_group_names)
  # Indexes of rows in cell_info_df that belong to duplicated groups
  inds_dupl_group_rows_lst <- lapply(dupl_group_names,
                                     function(x) which(cell_info_df$group == x))
  names(inds_dupl_group_rows_lst) <- dupl_group_names
  
  # Combination of representatives from each duplicated group
  represent_dupl_group_combis_df <- expand.grid(inds_dupl_group_rows_lst)
  return(
    list(represent_dupl_group_combis_df = represent_dupl_group_combis_df,
         unique_group_abbrs = unique_group_abbrs,
         dupl_group_names = dupl_group_names)
  )
  
}
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Get 1 representative per custom group of cell types
cell_info_df_orig <- read.csv(cell_info_path)
cell_info_df <- cell_info_df_orig

select_represent_group_out <- select_represent_group(cell_info_df)
represent_dupl_group_combis_df <- select_represent_group_out$represent_dupl_group_combis_df
unique_group_abbrs <- select_represent_group_out$unique_group_abbrs
dupl_group_names <- select_represent_group_out$dupl_group_names
combis_len <- nrow(represent_dupl_group_combis_df)

# List of groups of cell types to drop for each recalculation of Cp (next script)
drop_similar_lst <- lapply(1:combis_len, function(i) {
  
  represent_dupl_group_combis <- unname(unlist(represent_dupl_group_combis_df[i, ]))
  represent_dupl_group_abbrs <- cell_info_df$abbr[represent_dupl_group_combis]
  drop_abbrs <- setdiff(cell_info_df$abbr,
                        c(unique_group_abbrs, represent_dupl_group_abbrs))
  
  # Check that cell types (abbr) to be dropped belong to duplicated groups
  df_tmp <- cell_info_df[!cell_info_df$abbr %in% drop_abbrs, ]
  if (any(duplicated(df_tmp$group))) {
    stop("Combi", i, ": >1 dataset in a group.")
  }
  
  drop_groups <- cell_info_df$group[cell_info_df$abbr %in% drop_abbrs]
  if (!all(drop_groups %in% dupl_group_names)) {
    stop("Combi", i, ": Dropping cell types not in duplicated group.")
  } else {
    return(drop_abbrs)
  }
  
})

names(drop_similar_lst) <- paste0("drop_similar-", 1:length(drop_similar_lst))

# Identify datasets with outlying number of read pairs

dataqc_df <- read.csv(dataqc_path)
# Should be all reads per dataset and not just a certain scale because difference at a certain scale might be because of biological variation
read_type <- "mapped.monoclonal.read.pairs" #"mapped.monoclonal.read.pairs" #"cis" #"cis.15kb"

drop_read_outliers_lst <- list()

#bp_stats <- boxplot.stats(dataqc_df[[read_type]], coef = 1.5)$stats
#lower_bound <- bp_stats[1]
#upper_bound <- bp_stats[5]
#is_outlier <- dataqc_df[[read_type]] < lower_bound | dataqc_df[[read_type]] > upper_bound
#drop_read_outliers_lst[[paste0("drop_read_outliers-1.5iqr")]] <- dataqc_df$abbr[is_outlier]

for (n_mads in 3:1) {
  is_outlier <- scuttle::isOutlier(dataqc_df[[read_type]], nmads = n_mads, "both")
  drop_read_outliers_lst[[paste0("drop_read_outliers-", n_mads, "nmads")]] <- dataqc_df$abbr[is_outlier]
}

# Drop similar cell types and outlying number fo read pairs

lst_names <- setNames(c("drop_similar_outliers_2nmads-", "drop_similar_outliers_1nmads-"),
                      nm = c("drop_read_outliers-2nmads", "drop_read_outliers-1nmads"))

drop_similar_outliers_lst <- list(list(), list())
names(drop_similar_outliers_lst) <- unname(lst_names)

for (outlier_id in names(lst_names)) {
  
  read_outliers_abbs = drop_read_outliers_lst[[outlier_id]]
  cell_info_df <- cell_info_df_orig[!cell_info_df_orig$abbr %in% read_outliers_abbs,]
  
  represent_dupl_group_combis_df <- select_represent_group(cell_info_df)$represent_dupl_group_combis_df
  combis_len <- nrow(represent_dupl_group_combis_df)
  
  # List of groups of cell types to drop for each recalculation of Cp (next script)
  nme <- lst_names[outlier_id]
  drop_similar_outliers_lst[[nme]] <- lapply(1:combis_len, function(i) {
    
    c(cell_info_df$abbr[as.numeric(represent_dupl_group_combis_df[i,])], 
      read_outliers_abbs)
    
  })
  
  # Check
  len <- length(drop_similar_outliers_lst[[nme]])
  for (i in 1:len) {
    
    df <- cell_info_df_orig[!cell_info_df_orig$abbr %in% drop_similar_outliers_lst[[nme]][[i]],]
    if (any(duplicated(df$group) | df$abbr %in% read_outliers_abbs)) {
      stop("Error in ", outlier_id, " element ", i)
    }
          
  }
  
}

drop_similar_outliers_lst <- unlist(drop_similar_outliers_lst, recursive = FALSE)

# Assemble leave_out_lst

drop_one_lst <- setNames(as.list(celltype_ids), paste0("drop_one-", celltype_ids))
leave_out_lst <- c(drop_one_lst,
                   drop_similar_lst,
                   drop_read_outliers_lst,
                   drop_similar_outliers_lst)

leave_out_sorted_lst <- lapply(leave_out_lst, sort)
if (any(duplicated(leave_out_sorted_lst))) {
  stop("Duplicated element in leave_out_lst")
}

# REMOVE
# set.seed(290)
# leave_out_lst <- c(setNames(as.list(celltype_ids), paste0("drop_one-", celltype_ids))[1:2],
#                    drop_similar_lst[1:2],
#                    drop_read_outliers_lst,
#                    drop_similar_outliers_lst[sample(1:length(drop_similar_outliers_lst),size = 5)])

saveRDS(leave_out_lst, 
        file.path(wk.dir, "out_leave_out_list", "leave_out_list.rds"))

# rm(list=ls()); gc()
