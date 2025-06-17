# Set environment
library(tidyverse)
setwd("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow")
source("scripts/workflow/templates/canme_strains.r")

# Load BED files of LOH blocks with the allele
loh_alleles_files <- list.files(
    "results/annot_loh",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "LOH_regions.bed$"
)

names(loh_alleles_files) <- gsub(".*/", "", loh_alleles_files)
names(loh_alleles_files) <- gsub("\\..*", "", names(loh_alleles_files))
loh_alleles_files <- loh_alleles_files[order_eStrains_canme]

loh_alleles_list <- sapply(loh_alleles_files, FUN = function(file){
    if (file.size(file) == 0L) { 
        print("Empty file")
    } else {
        read.table(file, sep = "\t", header = FALSE)
    }
}, simplify = FALSE)

# Load GFF of genes that intersect at any point with an LOH block
# larger than 100bp and not in T0
loh_genes_files <- list.files(
    "results/annot_loh", 
    recursive = TRUE,
    full.names = TRUE,
    pattern = "LOH_genes.gff$"
)

names(loh_genes_files) <- gsub(".*/", "", loh_genes_files)
names(loh_genes_files) <- gsub("\\..*", "", names(loh_genes_files))
loh_genes_files <- loh_genes_files[order_eStrains_canme]

loh_genes_list <- sapply(loh_genes_files, FUN = function(file){
    read.table(file, sep = "\t", header = FALSE, comment.char = "") %>%
        # Filter mRNAs because they have the Parent data
        filter(V3 == "mRNA") %>%
        mutate(cds_length = V5 - V4) %>%
        # Split misc data
        mutate(misc_list = strsplit(V9, ";")) %>%
        mutate(misc_split = lapply(misc_list, function(str){
            split <- strsplit(str, "=")
            names(split) <- sapply(split, "[", 1)
            split <- lapply(split, "[", 2)
        })) %>%
        # Extract relevant fields
        mutate(geneID = sapply(misc_split, "[[", "ID")) %>%
        mutate(Preferred_name = sapply(misc_split, "[[", "em_Preferred_name")) %>%
        replace(.=="NULL", NA) %>%
        mutate(Preferred_name = unlist(Preferred_name)) %>%
        mutate(Parent = sapply(misc_split, "[[", "Parent")) %>%
        replace(.=="NULL", NA) %>%
        mutate(Parent = unlist(Parent)) %>%
        select(V1, V4, V5, geneID, misc_split, Preferred_name, Parent, cds_length) %>%
        distinct()
}, simplify = FALSE)

# Load GFF of the exact regions that intersect with LOH blocks
# larger than 100bp and not in T0
loh_region_files <- list.files(
    "results/annot_loh", 
    recursive = TRUE,
    full.names = TRUE,
    pattern = "LOH_region.gff$"
)

names(loh_region_files) <- gsub(".*/", "", loh_region_files)
names(loh_region_files) <- gsub("\\..*", "", names(loh_region_files))
loh_region_files <- loh_region_files[order_eStrains_canme]

loh_region_list <- sapply(loh_region_files, FUN = function(file){
    read.table(file, sep = "\t", header = FALSE, comment.char = "") %>%
        filter(V3 == "mRNA") %>%
        mutate(region_length = V5 - V4) %>%
        mutate(misc_list = strsplit(V9, ";")) %>%
        mutate(misc_split = lapply(misc_list, function(str){
            split <- strsplit(str, "=")
            names(split) <- sapply(split, "[", 1)
            split <- lapply(split, "[", 2)
        })) %>%
        mutate(geneID = sapply(misc_split, "[[", "ID")) %>%
        mutate(Preferred_name = sapply(misc_split, "[[", "em_Preferred_name")) %>%
        replace(.=="NULL", NA) %>%
        mutate(Preferred_name = unlist(Preferred_name)) %>%
        mutate(Parent = sapply(misc_split, "[[", "Parent")) %>%
        replace(.=="NULL", NA) %>%
        mutate(Parent = unlist(Parent)) %>%
        select(V1, V4, V5, geneID, Preferred_name, Parent, region_length)
}, simplify = FALSE)

# Melt alleles
molten_loh_alleles <- loh_alleles_list %>%
    map_df(~ .x, .id = "eStrain") %>%
    rename(all_of(c(chromosome = "V1", start = "V2", end = "V3", allele = "V4"))) %>%
    mutate(start = start + 1) %>%
    unique()

# Melt LOH genes
molten_loh_genes <- loh_genes_list %>%
    map_df(~ .x, .id = "eStrain") %>%
    mutate(eStrain_fct = factor(eStrain,
        levels = order_eStrains_canme)) %>%
    mutate(cond = conditions_canme[eStrain_fct]) %>%
    mutate(strain = case_when(
        eStrain %in% BP57_strains ~ "BP57",
        eStrain %in% DMKU_13A_strains ~ "DMKU_13A")) %>%
    group_by(eStrain, geneID) %>%
    mutate(V4 = min(V4)) %>%
    mutate(V5 = max(V5)) %>%
    mutate(gene_length = V5 - V4) %>%
    ungroup() %>%
    select(-cds_length) %>%
    distinct()

# Melt LOH regions
molten_loh_regions <- loh_region_list %>%
    map_df(~ .x, .id = "eStrain")  %>%
    mutate(eStrain_fct = factor(eStrain,
        levels = order_eStrains_canme)) %>%
    mutate(cond = conditions_canme[eStrain_fct]) %>%
    mutate(strain = case_when(
        eStrain %in% BP57_strains ~ "BP57",
        eStrain %in% DMKU_13A_strains ~ "DMKU_13A")) %>%
    distinct()

# Join each region with its corresponding gene feature
molten_loh_all <- left_join(molten_loh_genes, molten_loh_regions,
    by = c("geneID", "eStrain"), suffix = c("", "")) %>%
    ungroup() %>%
    rename(all_of(c(chromosome = "V1", start = "V4", end = "V5"))) %>%
    mutate(Preferred_name = case_when(
    Preferred_name == "-" | is.na(Preferred_name) ~ Parent,
    .default = Preferred_name
    )) %>%
    select(-misc_split)

# Join each region with its corresponding allele
molten_loh_wallele <-  inner_join(molten_loh_alleles, molten_loh_all,
    by = c("eStrain", "chromosome", "start", "end"),
    suffix = c("", "")) %>%
    ungroup() %>%
    unique()

# Save
write.csv(
    file = "results/annot_plots/loh_data_canme.csv",
    molten_loh_wallele,
    row.names = FALSE
)