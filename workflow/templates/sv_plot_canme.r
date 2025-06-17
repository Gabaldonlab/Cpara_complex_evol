# Set environment
library(tidyverse)

# Source the script to set the environment
source("scripts/workflow/templates/canme_strains.r")

# Load data
snp_data_canme <- read.csv("results/annot_plots/svs_data_canme.csv")

cat("Preparing data\n")
# Generate the data for the plots
molten_data <- snp_data_canme %>%
    mutate(SVtype = sapply(strsplit(X.Uploaded_variation, "\\|"), "[[", 1)) %>%
    mutate(chromosome = sapply(strsplit(Location, ":"), "[[", 1)) %>%
    mutate(positions = sapply(strsplit(Location, ":"), "[[", 2)) %>%
    mutate(start_pos = as.integer(sapply(strsplit(positions, "-"), "[[", 1))) %>%
    mutate(end_pos = as.integer(sapply(positions, FUN = function(pos){
        if (grepl("-", pos)) {
            strsplit(pos, "-")[[1]][2]
        } else {
            pos
        }
    }))) %>%
    mutate(eStrain = factor(eStrain,
        levels = order_eStrains_canme)) %>%
    mutate(strain = case_when(
        eStrain %in% BP57_strains ~ "BP57",
        eStrain %in% DMKU_13A_strains ~ "DMKU_13A")) %>%
    # Remove BND from the consequences
    mutate(Consequence = gsub("_BND", "", Consequence)) %>%
    # Split the consequences
    mutate(Consequence = strsplit(Consequence, ",")) %>%
    unnest_longer(Consequence) %>%
    # Mark inherited/new/reverted SNPs in switch conditions
    rowwise(eStrain) %>%
    mutate(switch_group = case_when(
        eStrain %in% names(switch_eStrains) ~ switch_eStrains[levels(eStrain)[eStrain]],
        .default = eStrain
    )) %>%
    group_by(SVtype, em_Preferred_name, switch_group) %>%
    mutate(switch = case_when(
        eStrain %in% switch_eStrains & any(eStrain %in% names(switch_eStrains)) ~ "inherited",
        eStrain %in% names(switch_eStrains) & all(!(eStrain %in% switch_eStrains)) ~ "reverted",
        eStrain %in% switch_eStrains ~ "new",
        eStrain %in% names(switch_eStrains) ~ "new",
        .default = NULL
    )) %>%
    ungroup() %>%
    # Sort by start position of the feature
    arrange(Start) %>%
    mutate(em_Preferred_name = factor(em_Preferred_name,
        levels = unique(em_Preferred_name))) %>%
    # Filter features that start before the start breakpoint (intersect with it)
    # or that end after the end breakpoint
    filter(Start < start_pos | end_pos < End)


# Colours for plots
# For shapes, filled shapes for sure calls and empty shapes for -like calls
SVtpye_shapes = c(
    # These are the tags used by CLOVE, DEL pairs with TDUP
    "CVD" = 17, "DEL" = 3, "DELlike" = 3, "TANlike" = 2, "TDUP" = 17,
    # These are genuine coverage-based deletions and duplications
    "coverageDEL" = 4, "coverageDUP" = 17,
    "INS" = 25,
    # These are in-chromosome inversions
    "INV" = 18, "INV1like" = 5, "INV2like" = 5,
    # These are all translocations between chromosomes
    "TRA" = 15, "IVT" = 15,
    "INVTX1like" = 0, "INVTX2like" = 0, "ITX1like" = 0, "ITX2like" = 0
)

cat("Plotting\n")
mutated_gene_plot <- ggplot(
    molten_data,
    aes(x = eStrain,
        y = em_Preferred_name,
        fill = Consequence,
        shape = SVtype,
        colour = Consequence,
        alpha = switch)) +
    geom_jitter(height = 0) +
    facet_grid(chromosome ~ cond, scales = "free", space = "free") +
    scale_y_discrete(labels = unique(molten_data)) +
    scale_alpha_manual(
        values = c("new" = 1, "inherited" = 0.1, "reverted" = 0.5),
        na.translate = FALSE
    ) +
    scale_shape_manual(
        values = SVtpye_shapes
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 2),
        shape = guide_legend(ncol = 2),
        alpha = guide_legend(ncol = 1)) +
    labs(title = "Protein-altering structural variants in evolved hybrids",
        subtitle = "All annotated genes",
        x = "Strain", y = "Gene")

ggsave("results/annot_plots/canme_all_svs.pdf", mutated_gene_plot,
    width = 10, height = 50, limitsize = FALSE)

# Generate plot only with the selected genes
select_gene_plot <- ggplot(
    molten_data %>% 
        # Filter out unnamed genes
        filter(em_Preferred_name != "") %>%
        # Subset selected genes
        filter(toupper(em_Preferred_name) %in% selected_genes),
    aes(x = eStrain,
        y = snp,
        fill = Consequence,
        shape = SVtype,
        alpha = switch)) +
    geom_tile(colour = "gray70") +
    scale_alpha_manual(
        values = c("new" = 1, "inherited" = 0.1, "reverted" = 0.5),
        na.translate = FALSE
    ) +
    facet_grid(em_Preferred_name ~ cond, scales = "free", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 2)) +
    labs(title = "Protein-altering variants in evolved hybrids",
        subtitle = "Selected genes",
        x = "Strain", y = "Mutation")

ggsave("results/annot_plots/canme_selected_svs.pdf", select_gene_plot,
    width = 6, height = 7, limitsize = FALSE)