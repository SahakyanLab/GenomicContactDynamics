###############################################################################
# HiC21 cell/tissue codes
################################################################################
CTCODE.MX <- rbind(c("Lymphoblast", "GM12878", "LC"),
                 c("Human_Embryonic_Stem_Cell", "H1", "ESC"),
                 c("Fetal_Lung_Fibroblast", "IMR90", "FC"),
                 c("Fetal_Lung_Fibroblast", "IMR-90", "FC"),
                 c("Mesendoderm", "MES", "MesC"),
                 c("Mesenchymal_Stem_Cell", "MSC", "MSC"),
                 c("Neural_Progenitor_Cell", "NPC", "NPC"),
                 c("Trophoblast-like_Cell", "TRO", "TLC"),
                 c("Adrenal_Gland", "AD", "AG"),
                 c("Aorta", "AO", "Ao"),
                 c("Bladder", "BL", "Bl"),
                 c("Dorsolateral_Prefrontal_Cortex", "CO", "Co"),
                 c("Hippocampus", "HC", "Hi"),
                 c("Lung", "LG", "Lu"),
                 c("Liver", "LI", "Li"),
                 c("Left_Ventricle", "LV", "LV"),
                 c("Ovary", "OV", "Ov"),
                 c("Pancreas", "PA", "Pa"),
                 c("Psoas_Muscle", "PO", "PM"),
                 c("Right_Ventricle", "RV", "RV"),
                 c("Small_Bowel", "SB", "SB"),
                 c("Spleen", "SX", "Sp") )
dimnames(CTCODE.MX) <- list(NULL, c("Celltiss", "Orig.code", "New.code"))
