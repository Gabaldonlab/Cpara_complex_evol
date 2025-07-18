# Set environment
library(tidyverse)
setwd("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow")
source("scripts/workflow/templates/canor_strains.r")

# Load data
molten_svs_all <- read.csv("results/annot_plots/svs_data_canor.csv")

# Remove up/downstream consequences
molten_svs_all <- molten_svs_all %>%
    # Subset selected genes
    #filter(toupper(em_Preferred_name) %in% selected_genes) %>%
    filter(!(Consequence %in% c(
        "upstream_gene_variant_BND", "upstream_gene_variant",
        "downstream_gene_variant_BND", "downstream_gene_variant"
    ))) %>%
    separate_rows(Consequence, sep = ",") %>%
    mutate(Consequence = gsub("_BND", "", Consequence),
        Consequence_simple = case_when(
            Consequence == "3_prime_UTR_variant" ~ "UTR-Intron",
            Consequence == "5_prime_UTR_variant" ~ "UTR-Intron",
            Consequence == "intron_variant" ~ "UTR-Intron",
            Consequence == "start_lost" ~ "Functional loss",
            Consequence == "splice_region_variant" ~ "Functional loss",
            Consequence == "coding_sequence_variant" ~ "CDS",
            Consequence == "protein_altering_variant" ~ "Missense",
            Consequence == "frameshift_variant" ~ "Missense",
            Consequence == "stop_gained" ~ "Missense",
            Consequence == "inframe_insertion" ~ "Insertion",
            Consequence == "feature_elongation" ~ "Insertion",
            Consequence == "transcript_amplification" ~ "Insertion",
            Consequence == "inframe_deletion" ~ "Deletion",
            Consequence == "feature_truncation" ~ "Deletion",
            Consequence == "transcript_ablation" ~ "Deletion",
            .default = Consequence
        )) %>%
    unique()

# Split 

brkpt_gene_plot <- ggplot(
    molten_svs_all,
    aes(x = eStrain,
        y = em_Preferred_name,
        fill = SVTYPE)) +
    geom_tile(colour = "gray70", alpha = 0.3) +
    geom_jitter(aes(colour = Consequence)) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 1), colour = guide_legend(ncol = 1)) +
    labs(title = "SV breakpoints affecting genes in evolved hybrids",
        subtitle = "Interest genes",
        x = "Strain", y = "Gene") +
    facet_grid(Scaffold ~ ., scale = "free_y", space = "free_y")

ggsave("results/annot_plots/canor_all_breakpoints.pdf", brkpt_gene_plot,
    width = 20, height = 55, limitsize = FALSE)