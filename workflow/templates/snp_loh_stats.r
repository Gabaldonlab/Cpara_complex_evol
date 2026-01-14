# Read files
# Set environment
library(tidyverse)
library(data.table)
setwd("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow")

# Select annotation tables and vcf files
lohalt_files <- list.files(
    "results/variant_loh",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "lohAltVars"
)
lohref_files <- list.files(
    "results/variant_loh",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "lohRefVars"        
)
nonloh_files <- list.files(
    "results/variant_loh",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "nonLohVars"        
)

parental_strains <- c(
    "EPOM_T0_7", "EPOM_T0_12", "EPOM_075_A7F", "EPOM_267_AinF7F",
    "EPOM_093_A12A", "EPOM_094_A12C", "EPOM_095_A12E", "EPOM_141_F12A",
    "EPOM_144_F12G", "EPOM_237_FinA12A", "EPOM_285_AinF12A"
)
# Select files
lohalt_match <- unlist(sapply(
    parental_strains,
    FUN = function(name){grep(name, lohalt_files)}
))
lohref_match <- unlist(sapply(
    parental_strains,
    FUN = function(name){grep(name, lohref_files)}
))
nonloh_match <- unlist(sapply(
    parental_strains,
    FUN = function(name){grep(name, nonloh_files)}
))
lohalt_files <- lohalt_files[lohalt_match]
lohref_files <- lohref_files[lohref_match]
nonloh_files <- nonloh_files[nonloh_match]

# Read and parse names
lohalt_list <- sapply(lohalt_files, FUN = function(file){
    fread(file, sep = "\t", header = TRUE, col.names =
        c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")) %>%
        rowwise() %>%
        mutate(GT = strsplit(SAMPLE, ":")[[1]][1])
}, simplify = FALSE)
names(lohalt_list) <- gsub("_lohAltVars.*", "", gsub(".*/", "", names(lohalt_list)))
lohalt_data <- lohalt_list %>% bind_rows(.id = "eStrain") %>% mutate(loh = "lohAlt")

lohref_list <- sapply(lohref_files, FUN = function(file){
    fread(file, sep = "\t", header = TRUE, col.names =
        c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")) %>%
        rowwise() %>%
        mutate(GT = strsplit(SAMPLE, ":")[[1]][1])
}, simplify = FALSE)
names(lohref_list) <- gsub("_lohRefVars.*", "", gsub(".*/", "", names(lohref_list)))
lohref_data <- lohref_list %>% bind_rows(.id = "eStrain") %>% mutate(loh = "lohRef")

nonloh_list <- sapply(nonloh_files, FUN = function(file){
    fread(file, sep = "\t", header = TRUE, col.names =
        c("CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")) %>%
        rowwise() %>%
        mutate(GT = strsplit(SAMPLE, ":")[[1]][1])
}, simplify = FALSE)
names(nonloh_list) <- gsub("_nonLohVars.*", "", gsub(".*/", "", names(nonloh_list)))
nonloh_data <- nonloh_list %>% bind_rows(.id = "eStrain") %>% mutate(loh = "nonLoh")

allsnp_data <- rbind(lohalt_data, lohref_data, nonloh_data)

lohratio_per_chr <- allsnp_data %>% 
    group_by(eStrain, loh, CHROM) %>%
    summarise(n = n()) %>%
    pivot_wider(id_cols = c(eStrain, CHROM), names_from = loh, values_from = n) %>%
    mutate(lohAlt = case_when(is.na(lohAlt) ~ 0, .default = lohAlt)) %>%
    mutate(lohRef = case_when(is.na(lohRef) ~ 0, .default = lohRef)) %>%
    mutate(nonLoh = case_when(is.na(nonLoh) ~ 0, .default = nonLoh)) %>%
    mutate(loh_ratio = (lohAlt + lohRef) / (lohAlt + lohRef + nonLoh) ) %>%
    arrange(desc(loh_ratio)) 
    
lohratio_global <- lohratio_per_chr %>%
    group_by(eStrain) %>%
    summarise(
        lohAlt = sum(lohAlt),
        lohRef = sum(lohRef),
        nonLoh = sum(nonLoh)) %>%
    mutate(loh_ratio = (lohAlt + lohRef) / (lohAlt + lohRef + nonLoh) ) %>%
    arrange(desc(loh_ratio)) 


    %>%
    pivot_wider(id_cols = c(eStrain), names_from = CHROM, values_from = loh_ratio)
