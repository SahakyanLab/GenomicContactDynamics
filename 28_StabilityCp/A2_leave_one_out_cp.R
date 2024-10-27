################################################################################
# Recalculate Cp after removing a cell type or group of cell types, scale Cp by
# number of celltype/tissue used
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
data.dir = file.path(home.dir, "data/gcd/Database")
persist.dir = file.path(data.dir, "HiC_features_GSE87112_RAWpc")
leave_out_list_path = file.path(wk.dir, "out_leave_out_list/leave_out_list.rds")
out.dir = paste0(wk.dir, "/out_leave_one_out_cp")
### OTHER SETTINGS #############################################################
celltype_ids_ref = c("Co", "Hi", "Lu", "LV", "RV", "Ao", "PM", "Pa", "Sp", "Li", "SB",
                     "AG", "Ov", "Bl", "MesC", "MSC", "NPC", "TLC", "ESC", "FC", "LC")
gcb = "min2Mb"
chrs = paste0("chr", c("X", 1:22))
subset_leave_out = NULL #c("drop_abbrs_outliers_1.5iqr", "drop_abbrs_outliers_3nmads", "drop_abbrs_outliers_2nmads")
save_orig_cp = TRUE #FALSE
################################################################################
# LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES * LIBRARIES & DEPENDENCIES 
################################################################################
################################################################################
# MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE * MAIN CODE *
################################################################################
leave_out_lst <- readRDS(leave_out_list_path)
if (!is.null(subset_leave_out)) {
  leave_out_lst <- leave_out_lst[subset_leave_out]
}

leave_out_lst_len <- length(leave_out_lst)
if (any(duplicated(names(leave_out_lst)) | duplicated(leave_out_lst))) {
 stop("Duplicated elements in leave_out_lst") 
}

for (chr in chrs) {
  
  load(file.path(persist.dir, paste0(chr, "_Persist_", gcb, ".RData")))
  celltype_ids <- setdiff(colnames(PERSIST.MX$hits), c("i", "j"))
  
  if (!identical(celltype_ids_ref, celltype_ids)) {
    stop(gcb, " ", chr, ": celltype_ids not matching with celltype_ids_ref")
  }
  
  #recalc_cp_mx <- matrix(data = NA, nrow = nrow(PERSIST.MX$hits), ncol = leave_out_lst_len,
  #                       dimnames = list(rownames(PERSIST.MX$hits), names(leave_out_lst)))
  
  for (drop_id in names(leave_out_lst)) {
    
    drop_celltypes <- leave_out_lst[[drop_id]]
    mx <- data.matrix(PERSIST.MX$hits[, setdiff(celltype_ids, drop_celltypes), drop = FALSE])
    #recalc_cp_mx[,i] <- unname(rowSums(mx > 0)) / ncol(mx)
    
    recalc_cp <- unname(rowSums(mx > 0)) / ncol(mx)
    saveRDS(recalc_cp, file.path(out.dir,
                                 paste(drop_id, gcb, chr, "recalculated_cp.rds", sep = "_")))
    rm(mx, recalc_cp, drop_celltypes)
  }
  
  if (save_orig_cp) {
    
    orig_cp <- PERSIST.MX$ntis / length(celltype_ids)
    saveRDS(orig_cp, file.path(out.dir, paste(gcb, chr, "original_cp.rds", sep = "_")))
    rm(orig_cp)
    
  }
  
  #recalc_cp_mx <- cbind(recalc_cp_mx, orig_cp = PERSIST.MX$ntis / length(celltype_ids))
  #saveRDS(recalc_cp_mx, file.path(out.dir, paste0(gcb, "_", chr, "_recalculated_cp.rds")))
  
  message(gcb, " ", chr, ": done...")
  rm(PERSIST.MX, celltype_ids)
  
}

# rm(list=ls()); gc()

## Use package with p-value
#cor_mx <- cor(recalc_cp_mx)
#range(cor_mx[,"orig_cp"])

