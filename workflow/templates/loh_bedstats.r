# Libraries
library(ggpubr)
library(tidyverse)
library(data.table)

# Load metadata
metadata <- read.table("data/metadata/estrain_metadata_workflow.tsv",
    sep = "\t",
    header = TRUE)

# Load the Bedtools stats
file_vec <- list.files(
    "results/jloh_bedstats", 
    recursive = TRUE,
    full.names = TRUE,
    pattern = "_stats.tsv")

bedstats_list <- sapply(file_vec, FUN = function(file){
    fread(file, col.names = c(
        "chrom", "chrom_length", "num_ivls", "total_ivl_bp",
        "chrom_frac_genome", "frac_all_ivls", "frac_all_bp",
        "min", "max", "mean", "NA"
    ), fill = TRUE)
}, simplify = FALSE)

names(bedstats_list) <- sapply(file_vec,
    FUN = function(name){
        paste(strsplit(gsub(".*/", "", name), "_")[[1]][c(1,2,3)], collapse = "_")
    })

bedstats_df <- lapply(bedstats_list, FUN = function(tab){
    tab <- tab %>%
        mutate(ivls_frac_chrom = total_ivl_bp / chrom_length)
    }) %>%
    bind_rows(.id = "full_name") %>%
    filter(!(chrom %in% c("DQ026513.2", "JQ062879.1"))) %>%
    filter(full_name != "EPOM_224_FinA8G") %>%
    # Filter out the strains with whole genome LOH
    filter(!(full_name %in% c("EPOM_088_A10G", "EPOM_133_F10A", "EPOM_183_AF10E",
        "EPOM_229_FinA10A", "EPOM_280_AinF10G"))) %>%
    left_join(metadata) %>%
    mutate(evolved = case_when(
        condition == "Y" ~ FALSE,
        condition != "Y" ~ TRUE
    )) %>%
    filter(!(parental %in% c("Co90125", "SY36"))) %>%
    group_by(parental) %>%
    mutate(ypd_ivls_frac = ivls_frac_chrom[condition == "Y" & chrom == "all"] ) %>%
    group_by(full_name) %>%
    mutate(ivls_above_ypd = ivls_frac_chrom[chrom == "all"] > 0.005)

bedstats_all <- bedstats_df %>% 
    filter(chrom == "all")

# Plot the statistics
plot <- ggplot(
    bedstats_df %>% filter(chrom != "all"),
    aes(x = full_name, y = ivls_frac_chrom)) +
    #geom_point(aes(colour = ivls_above_ypd)) +
    geom_point(data = bedstats_all, colour = "firebrick", size = 2, shape = 1) +
    geom_hline(
        data = bedstats_df %>% filter(chrom == "all", condition == "Y"),
        aes(yintercept = ivls_frac_chrom)
    ) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    facet_grid(. ~ parental, space = "free", scales = "free")

ggsave("results/jloh_plot/loh_bedstats.pdf", plot, width = 9)

# Get statistics

# Average LOH in evolved strains
bedstats_all %>% filter(condition != "Y") %>% pull(ivls_frac_chrom) %>% mean()
bedstats_all %>% filter(condition != "Y") %>% pull(ivls_frac_chrom) %>% sd()
bedstats_all %>% filter(condition != "Y") %>% nrow()
# Average LOH in YPD strains
bedstats_all %>% filter(condition == "Y") %>% pull(ivls_frac_chrom) %>% mean()
bedstats_all %>% filter(condition == "Y") %>% pull(ivls_frac_chrom) %>% sd()
bedstats_all %>% filter(condition == "Y") %>% nrow()
# Average LOH in per speices
bedstats_all %>% filter(condition != "Y", reference == "CANOR_CO90125_mit_MCO471") %>% pull(ivls_frac_chrom) %>% mean()
bedstats_all %>% filter(condition != "Y", reference == "CANOR_CO90125_mit_MCO471") %>% pull(ivls_frac_chrom) %>% sd()
bedstats_all %>% filter(condition != "Y", reference == "CANOR_CO90125_mit_MCO471") %>% nrow()
bedstats_all %>% filter(condition != "Y", reference == "CANME_BP57_mit_BP57") %>% pull(ivls_frac_chrom) %>% mean()
bedstats_all %>% filter(condition != "Y", reference == "CANME_BP57_mit_BP57") %>% pull(ivls_frac_chrom) %>% sd()
bedstats_all %>% filter(condition != "Y", reference == "CANME_BP57_mit_BP57") %>% nrow()
# Test
wilcox.test(ivls_frac_chrom ~ evolved, bedstats_all)
wilcox.test(ivls_frac_chrom ~ reference, bedstats_all %>% filter(evolved))
