# Set environment
library(tidyverse)
library(data.table)
setwd("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow")
source("scripts/workflow/templates/canor_strains.r")
source("scripts/workflow/templates/canme_strains.r")

# Load the CNV regions
# This is the just the CNV region over the feature
# Contains the number of copies
file_vec <- list.files(
    "results/annot_cnv", 
    recursive = TRUE,
    full.names = TRUE,
    pattern = "CNV_GFF_merge.bed")

cnv_gff_intersect <- sapply(file_vec, FUN = function(file){
    print(file)
    if (file.size(file) != 0){
    read.table(file, sep = "\t", header = FALSE, comment.char = "",
        col.names = c("cnv_chr", "cnv_start", "cnv_end", "copy_number",
            "gene_chr", "source", "type", "gene_start", "gene_end", "score",
            "strand", "null", "misc", "overlap", "eStrain")) %>%
        # Caculate overlap
        mutate(species = case_when(
            eStrain %in% canor_strains ~ "CANOR",
            eStrain %in% order_eStrains_canme ~ "CANME"
        )) %>%
        mutate(gene_length = gene_end - gene_start - 1) %>%
        mutate(cnv_length = cnv_end - cnv_start) %>%
        mutate(perc_overlap = overlap / gene_length * 100) %>%
        mutate(cnv_type = case_when(
            copy_number > 1 ~ "DUP",
            copy_number < 1 ~ "DEL"
        )) %>%
        mutate(cnv_smaller_than_gene = cnv_length < gene_length) %>%
        # Split the miscellaneous column
        mutate(misc_list = strsplit(misc, ";")) %>%
        mutate(misc_split = lapply(misc_list, function(str){
            split <- strsplit(str, "=")
            names(split) <- sapply(split, "[", 1)
            split <- lapply(split, "[", 2)
        })) %>%
        mutate(geneID = sapply(misc_split, "[[", "ID")) %>%
        mutate(em_Preferred_name = sapply(misc_split, "[[", "em_Preferred_name")) %>%
        replace(.=="NULL", NA) %>%
        mutate(em_Preferred_name = unlist(em_Preferred_name)) %>%
        mutate(em_target = sapply(misc_split, "[[", "em_target")) %>%
        replace(.=="NULL", NA) %>%
        mutate(em_target = unlist(em_target)) %>%
        select(-c(misc_list, misc_split)) %>%
        mutate(gene_name = case_when(
            is.na(em_Preferred_name) ~ geneID,
            .default = em_Preferred_name
        ))
    }
}, simplify = FALSE)

molten_cnv_intersect <- cnv_gff_intersect %>%
    map_df(~.x)

filtered_cnv <- molten_cnv_intersect %>%
    # Remove duplications that are smaller than the gene
    filter(!cnv_smaller_than_gene | cnv_type == "DEL") %>%
    # Remove duplications that are not affected in more than 80%
    filter(perc_overlap > 80 | cnv_type == "DEL") %>%
    # Remove deletions shorter than 1kb
    filter(cnv_length > 1000 | cnv_type == "DUP") %>%
    # Remove deletions that are not affected in more than 80%
    filter(perc_overlap > 10 | cnv_type == "DUP") %>%
    # Filter the EPOM_224_FinA8G species
    filter(eStrain != "EPOM_224_FinA8G")

# Stats
# Number of unique "events"
filtered_cnv %>% 
    select(species, cnv_type, eStrain, cnv_chr, cnv_start, cnv_end, cnv_length) %>%
    # Keep one entry for each CNV, they have one for each affected gene
    unique() %>%
    # Summarise the number of events
    group_by(species, cnv_type) %>%
    summarise(
        n = n(),
        mean_length = mean(cnv_length),
        median_length = median(cnv_length))

# Number of genes affected
filtered_cnv %>% 
    group_by(species, cnv_type) %>%
    summarise(
        n_affected_genes = length(unique(geneID)),
        n_aff_annot_genes = length(unique(em_Preferred_name)))

# Wirte in the drive's format
write.table(filtered_cnv %>%
    arrange(species) %>%
    select(
        eStrain, gene_chr, gene_name, gene_start, gene_end, gene_length,
        misc, cnv_type, cnv_start, cnv_end, cnv_length, copy_number,
        overlap, perc_overlap, cnv_smaller_than_gene
    ),
    "results/annot_plots/cnv_data.tsv",
    sep = "\t",
    quote = FALSE,
    row.names = FALSE)