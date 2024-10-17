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
### OTHER SETTINGS #############################################################
cell_info_path = file.path(wk.dir, "celltypes.csv")
celltype_ids = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
                 "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
# Get 1 representative per custom group of cell types
cell_info_df <- read.csv(cell_info_path)

# Groups with only one cell type (abbr column in cell_info_df)
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
combis_len <- nrow(represent_dupl_group_combis_df)

# List of groups of cell types to drop for each recalculation of Cp (next script)
drop_abbrs_lst <- lapply(1:combis_len, function(i) {
  
  represent_dupl_group_combis <- unname(unlist(represent_dupl_group_combis_df[i, ]))
  represent_dupl_group_abbrs <- cell_info_df$abbr[represent_dupl_group_combis]
  drop_abbrs <- setdiff(cell_info_df$abbr,
                        c(unique_group_abbrs, represent_dupl_group_abbrs))
  
  # Check that cell types (abbr) to be dropped belong to duplicated groups
  drop_groups <- cell_info_df$group[cell_info_df$abbr %in% drop_abbrs]
  if (!all(drop_groups %in% dupl_group_names)) {
    stop("Combi", i, ": Dropping cell types not in duplicated group.")
  } else {
    return(drop_abbrs)
  }
  
})

names(drop_abbrs_lst) <- paste0("drop_abbrs_", 1:length(drop_abbrs_lst))
leave_out_lst <- c(setNames(as.list(celltype_ids), paste0("drop_", celltype_ids)), 
                   drop_abbrs_lst)
saveRDS(leave_out_lst, 
        file.path(wk.dir, "out_leave_out_list", "leave_out_list.rds"))
