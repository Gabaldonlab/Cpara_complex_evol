# Set env
BP57_strains <- c("EPOM_017_Y5B", "EPOM_066_A5D", "EPOM_068_A5H",
    "EPOM_113_F5B", "EPOM_115_F5F", "EPOM_116_F5H", "EPOM_209_FinA5B")

DMKU_13A_strains <- c("EPOM_021_Y6A", "EPOM_070_A6C", "EPOM_072_A6G",
    "EPOM_117_F6A", "EPOM_118_F6C", "EPOM_119_F6E", "EPOM_213_FinA6A",
    "EPOM_262_AinF6C", "EPOM_264_AinF6G")

order_eStrains_canme <- c(
            "EPOM_066_A5D", "EPOM_113_F5B", "EPOM_115_F5F", "EPOM_068_A5H", 
            "EPOM_116_F5H", "EPOM_209_FinA5B",
            "EPOM_117_F6A", "EPOM_118_F6C", "EPOM_119_F6E", 
            "EPOM_213_FinA6A", "EPOM_070_A6C", "EPOM_262_AinF6C",
            "EPOM_072_A6G", "EPOM_264_AinF6G")

order_allstrains_canme <- c(
            "EPOM_017_Y5B",
            "EPOM_066_A5D", "EPOM_068_A5H", "EPOM_115_F5F",
            "EPOM_116_F5H", "EPOM_113_F5B", "EPOM_209_FinA5B",
            "EPOM_021_Y6A",
            "EPOM_118_F6C", "EPOM_119_F6E", "EPOM_117_F6A",
            "EPOM_213_FinA6A", "EPOM_070_A6C", "EPOM_262_AinF6C",
            "EPOM_072_A6G", "EPOM_264_AinF6G")

ypd_strains_canme <- c(
    "BP57" = "EPOM_017_Y5B",
    "13A"  = "EPOM_021_Y6A"
)

conditions_canme <- c(
            "EPOM_017_Y5B" = "YPD", "EPOM_021_Y6A" = "YPD",
            "EPOM_066_A5D" = "ANI", "EPOM_068_A5H" = "ANI", "EPOM_115_F5F" = "FLZ",
            "EPOM_116_F5H" = "FLZ", "EPOM_113_F5B" = "FLZ", "EPOM_209_FinA5B" = "FinA",
            "EPOM_118_F6C" = "FLZ", "EPOM_119_F6E" = "FLZ", "EPOM_117_F6A" = "FLZ",
            "EPOM_213_FinA6A" = "FinA", "EPOM_070_A6C" = "ANI", "EPOM_262_AinF6C" = "AinF",
            "EPOM_072_A6G" = "ANI", "EPOM_264_AinF6G" = "AinF")

selected_genes <- c("UPC2", "UME6", "TAC1", "PDR16", "PDR1", "MRR1", "MDR1",
    "GSL2", "GSL1", "GSC1", "FKS1", "FKS2", "ERG4", "ERG3",
    "ERG2", "ERG13", "ERG12", "ERG11", "ERG1", "EFG1", "CZF1", "CPH2", "CHS3",
    "CDR1B", "CDR1", "BCR1")

switch_eStrains <- c(
    "EPOM_113_F5B" = "EPOM_209_FinA5B",
    "EPOM_117_F6A" = "EPOM_213_FinA6A",
    "EPOM_070_A6C" = "EPOM_262_AinF6C",
    "EPOM_072_A6G" = "EPOM_264_AinF6G"
)

genes_with_snps <- c(
    "MRR1",
    "FKS1",
    "ERG11",
    "NIK1",
    "GRE3",
    "rga1",
    "TOS9",
    "TGL5",
    "GLN4",
    "mRNA−g4364",
    "mRNA−g3224",
    "mRNA−g1798",
    "mRNA−g1480",
    "SRT1",
    "SCS7",
    "SAP98",
    "RRP5",
    "RHO1",
    "RFG1",
    "GPA2",
    "CCT5",
    "BUD5"
)

interest_genes <- unique(
    c(
        genes_with_snps,
        selected_genes
    )
)