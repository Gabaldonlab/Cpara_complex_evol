# Set env

# Parental A strains - CO90125 strain
pa_strains <- c(
    "EPOM_025_Y7B",
    "EPOM_075_A7F",
    "EPOM_267_AinF7F"
)

# Clade 1 strains - MCO456 strain
c1_strains <- c(
    "EPOM_029_Y8A",
    #"EPOM_224_FinA8G", excluded for not having switch
    "EPOM_079_A8E",
    "EPOM_271_AinF8E",
    "EPOM_078_A8C",
    "EPOM_270_AinF8C"
)

# Clade 2 strains - IFM48386 strain
c2_strains <- c(
    "EPOM_036_Y9H",
    "EPOM_129_F9B",
    "EPOM_225_FinA9B",
    "EPOM_132_F9H",
    "EPOM_228_FinA9H",
    "EPOM_083_A9F",
    "EPOM_275_AinF9F",
    "EPOM_084_A9H",
    "EPOM_276_AinF9H",
    "EPOM_177_AF9B"
)

# Clade 3 strains - CBS172 strain
c3_strains <- c(
    "EPOM_037_Y10A",
    "EPOM_134_F10C",
    "EPOM_230_FinA10C",
    "EPOM_133_F10A",
    "EPOM_229_FinA10A",
    "EPOM_088_A10G",
    "EPOM_280_AinF10G",
    "EPOM_087_A10E",
    "EPOM_183_AF10E"
)

# Clade 4 strains - MCO457 strain
c4_strains <- c(
    "EPOM_041_Y11B",
    "EPOM_138_F11D",
    "EPOM_234_FinA11D",
    "EPOM_139_F11F",
    "EPOM_235_FinA11F",
    "EPOM_089_A11B",
    "EPOM_281_AinF11B",
    "EPOM_091_A11F",
    "EPOM_283_AinF11F",
    "EPOM_186_AF11D"
)

# Parental B strains - SY36 strain
pb_strains <- c(
    "EPOM_045_Y12A",
    "EPOM_093_A12A",
    "EPOM_285_AinF12A",
    "EPOM_141_F12A",
    "EPOM_237_FinA12A",
    "EPOM_144_F12G",
    "EPOM_094_A12C",
    "EPOM_095_A12E"
)

# All strains
canor_strains <- c(
    pa_strains,
    c1_strains,
    c2_strains,
    c3_strains,
    c4_strains,
    pb_strains
)

# T0 strains 
t0_strains <- c(
    "pa_strains" = "EPOM_T0_7",
    "c1_strains" = "EPOM_T0_8",
    "c2_strains" = "EPOM_T0_9",
    "c3_strains" = "EPOM_T0_10",
    "c4_strains" = "EPOM_T0_11",
    "pb_strains" = "EPOM_T0_12"
)

# YPD strains
ypd_strains <- c(
    "pa_strains" = "EPOM_025_Y7B",
    "c1_strains" = "EPOM_029_Y8A",
    "c2_strains" = "EPOM_036_Y9H",
    "c3_strains" = "EPOM_037_Y10A",
    "c4_strains" = "EPOM_041_Y11B",
    "pb_strains" = "EPOM_045_Y12A"
)

# Suspicious haploids
hap_strains <- c(
    "EPOM_133_F10A",
    "EPOM_088_A10G",
    "EPOM_183_AF10E"
)

# Paired samples in switch conditions
switch_eStrains <- c(
    "EPOM_075_A7F" = "EPOM_267_AinF7F",
    "EPOM_079_A8E" = "EPOM_271_AinF8E",
    "EPOM_078_A8C" = "EPOM_270_AinF8C",
    "EPOM_129_F9B" = "EPOM_225_FinA9B",
    "EPOM_132_F9H" = "EPOM_228_FinA9H",
    "EPOM_083_A9F" = "EPOM_275_AinF9F",
    "EPOM_084_A9H" = "EPOM_276_AinF9H",
    "EPOM_134_F10C" = "EPOM_230_FinA10C",
    "EPOM_133_F10A" = "EPOM_229_FinA10A",
    "EPOM_088_A10G" = "EPOM_280_AinF10G",
    "EPOM_138_F11D" = "EPOM_234_FinA11D",
    "EPOM_139_F11F" = "EPOM_235_FinA11F",
    "EPOM_089_A11B" = "EPOM_281_AinF11B",
    "EPOM_091_A11F" = "EPOM_283_AinF11F",
    "EPOM_093_A12A" = "EPOM_285_AinF12A",
    "EPOM_141_F12A" = "EPOM_237_FinA12A"
)

# Order considering clade and evolution step
order_eStrain_canor <- order_eStrain <- c(
    "EPOM_075_A7F",
    "EPOM_267_AinF7F",
    #"EPOM_224_FinA8G", excluded for not having switch
    "EPOM_079_A8E",
    "EPOM_271_AinF8E",
    "EPOM_078_A8C",
    "EPOM_270_AinF8C",
    "EPOM_129_F9B",
    "EPOM_225_FinA9B",
    "EPOM_132_F9H",
    "EPOM_228_FinA9H",
    "EPOM_083_A9F",
    "EPOM_275_AinF9F",
    "EPOM_084_A9H",
    "EPOM_276_AinF9H",
    "EPOM_177_AF9B",
    "EPOM_134_F10C",
    "EPOM_230_FinA10C",
    "EPOM_133_F10A",
    "EPOM_229_FinA10A",
    "EPOM_088_A10G",
    "EPOM_280_AinF10G",
    "EPOM_087_A10E",
    "EPOM_183_AF10E",
    "EPOM_138_F11D",
    "EPOM_234_FinA11D",
    "EPOM_139_F11F",
    "EPOM_235_FinA11F",
    "EPOM_089_A11B",
    "EPOM_281_AinF11B",
    "EPOM_091_A11F",
    "EPOM_283_AinF11F",
    "EPOM_186_AF11D",
    "EPOM_093_A12A",
    "EPOM_285_AinF12A",
    "EPOM_141_F12A",
    "EPOM_237_FinA12A",
    "EPOM_144_F12G",
    "EPOM_094_A12C",
    "EPOM_095_A12E")

# Including YPD
order_allstrains_canor <- c(
    "EPOM_025_Y7B",
    "EPOM_075_A7F",
    "EPOM_267_AinF7F",
    "EPOM_029_Y8A",
    #"EPOM_224_FinA8G", excluded for not having switch  
    "EPOM_079_A8E",
    "EPOM_271_AinF8E",
    "EPOM_078_A8C",
    "EPOM_270_AinF8C",
    "EPOM_036_Y9H",
    "EPOM_129_F9B",
    "EPOM_225_FinA9B",
    "EPOM_132_F9H",
    "EPOM_228_FinA9H",
    "EPOM_083_A9F",
    "EPOM_275_AinF9F",
    "EPOM_084_A9H",
    "EPOM_276_AinF9H",
    "EPOM_177_AF9B",
    "EPOM_037_Y10A",
    "EPOM_134_F10C",
    "EPOM_230_FinA10C",
    "EPOM_133_F10A",
    "EPOM_229_FinA10A",
    "EPOM_088_A10G",
    "EPOM_280_AinF10G",
    "EPOM_087_A10E",
    "EPOM_183_AF10E",
    "EPOM_041_Y11B",
    "EPOM_138_F11D",
    "EPOM_234_FinA11D",
    "EPOM_139_F11F",
    "EPOM_235_FinA11F",
    "EPOM_089_A11B",
    "EPOM_281_AinF11B",
    "EPOM_091_A11F",
    "EPOM_283_AinF11F",
    "EPOM_186_AF11D",
    "EPOM_045_Y12A",
    "EPOM_093_A12A",
    "EPOM_285_AinF12A",
    "EPOM_141_F12A",
    "EPOM_237_FinA12A",
    "EPOM_144_F12G",
    "EPOM_094_A12C",
    "EPOM_095_A12E")

conditions_canor <- c(
    "EPOM_025_Y7B" = "YPD",
    "EPOM_029_Y8A" = "YPD",
    "EPOM_036_Y9H" = "YPD",
    "EPOM_037_Y10A" = "YPD",
    "EPOM_041_Y11B" = "YPD",
    "EPOM_045_Y12A" = "YPD",
    "EPOM_075_A7F" = "ANI",
    "EPOM_267_AinF7F" = "AinF",
    "EPOM_093_A12A" = "ANI",
    #"EPOM_224_FinA8G" = "FinA", excluded for not having switch
    "EPOM_079_A8E" = "ANI",
    "EPOM_271_AinF8E" = "AinF",
    "EPOM_078_A8C" = "ANI",
    "EPOM_270_AinF8C" = "AinF",
    "EPOM_129_F9B" = "FLZ",
    "EPOM_225_FinA9B" = "FinA",
    "EPOM_132_F9H" = "FLZ",
    "EPOM_228_FinA9H" = "FinA",
    "EPOM_083_A9F" = "ANI",
    "EPOM_275_AinF9F" = "AinF",
    "EPOM_084_A9H" = "ANI",
    "EPOM_276_AinF9H" = "AinF",
    "EPOM_177_AF9B" = "ANIFLZ",
    "EPOM_134_F10C" = "FLZ",
    "EPOM_230_FinA10C" = "FinA",
    "EPOM_133_F10A" = "FLZ",
    "EPOM_229_FinA10A" = "FinA",
    "EPOM_088_A10G" = "ANI",
    "EPOM_280_AinF10G" = "AinF",
    "EPOM_087_A10E" = "ANI",
    "EPOM_183_AF10E" = "ANIFLZ",
    "EPOM_138_F11D" = "FLZ",
    "EPOM_234_FinA11D" = "FinA",
    "EPOM_139_F11F" = "FLZ",
    "EPOM_235_FinA11F" = "FinA",
    "EPOM_089_A11B" = "ANI",
    "EPOM_281_AinF11B" = "AinF",
    "EPOM_091_A11F" = "ANI",
    "EPOM_283_AinF11F" = "AinF",
    "EPOM_186_AF11D" = "ANIFLZ",
    "EPOM_285_AinF12A" = "AinF",
    "EPOM_141_F12A" = "FLZ",
    "EPOM_237_FinA12A" = "FinA",
    "EPOM_144_F12G" = "FLZ",
    "EPOM_094_A12C" = "ANI",
    "EPOM_095_A12E" = "ANI")

selected_genes <- c(
    "UPC2",
    "UME6",
    "TAC1",
    "SSC1",
    "PDR16",
    "PDR1",
    "MRR1",
    "MDR1",
    "GSL2",
    "GSL1",
    "GSC1",
    "FKS1",
    "FKS2",
    "ERG4",
    "ERG3",
    "ERG2",
    "ERG13",
    "ERG12",
    "ERG11",
    "ERG1",
    "EFG1",
    "CZF1",
    "CPH2",
    "CHS3",
    "CDR1B",
    "CDR1",
    "BCR1"
)

genes_with_snps <- c(
    "BMT2",
    "CORT_0G03370",
    "GAL4",
    "YAP1",
    "CORT_0D00380",
    "HSX11",
    "FTH1",
    "SEN34",
    "RPL21A",
    "STE11",
    "KCH1",
    "CORT_0C01840",
    "HSR1",
    "NPR3",
    "CORT_0D04430",
    "CORT_0F00790",
    "CORT_0G02660",
    "TAC1",
    "ERG11",
    "SSC1",
    "FKS1",
    "MRR1"
)

interest_genes <- unique(
    c(
        genes_with_snps,
        selected_genes
    )
)