# Set environment
library(tidyverse)

# Source the script to set the environment
args <- commandArgs(trailingOnly = TRUE)
source(args[1])

# Load data
snp_data_canme <- read.csv("snp_data_canme.csv")

cat("Plotting\n")
# Generate a plot of mutations in genes
molten_data <- snp_data_canme %>%
    mutate(snp = paste0(Protein_position, Amino_acids)) %>%
    mutate(cond = conditions_canme[eStrain]) %>%
    mutate(eStrain = factor(eStrain,
        levels = order_eStrains_canme)) %>%
    mutate(strain = case_when(
        eStrain %in% BP57_strains ~ "BP57",
        eStrain %in% DMKU_13A_strains ~ "DMKU_13A")) %>%
    # Mark inherited/new/reverted SNPs in switch conditions
    ungroup() %>%
    rowwise(eStrain) %>%
    mutate(switch_group = case_when(
        eStrain %in% names(switch_eStrains) ~ switch_eStrains[levels(eStrain)[eStrain]],
        .default = eStrain
    )) %>%
    group_by(X.Uploaded_variation, switch_group) %>%
    mutate(switch = case_when(
        eStrain %in% switch_eStrains & any(eStrain %in% names(switch_eStrains)) ~ "inherited",
        eStrain %in% names(switch_eStrains) & all(!(eStrain %in% switch_eStrains)) ~ "reverted",
        eStrain %in% switch_eStrains ~ "new",
        eStrain %in% names(switch_eStrains) ~ "new",
        .default = NULL
    )) %>%
    # Get the number of mutations in each gene
    group_by(eStrain, Feature) %>%
    mutate(cumm_snps = length(unique(X.Uploaded_variation))) %>%
    group_by(Feature) %>%
    # Sort by number of strains with a mutation
    mutate(total_snps = length(unique(levels(eStrain)[eStrain]))) %>%
    arrange(total_snps)

mutated_gene_plot <- ggplot(
    molten_data,
    aes(x = eStrain,
        y = snp,
        fill = Consequence,
        alpha = switch)) +
    geom_tile(colour = "gray70") +
    facet_grid(em_Preferred_name ~ cond, scales = "free", space = "free") +
    scale_y_discrete(labels = unique(molten_data)) +
    scale_alpha_manual(
        values = c("new" = 1, "inherited" = 0.1, "reverted" = 0.5),
        na.translate = FALSE
    ) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 2)) +
    labs(title = "Protein-altering variants in evolved hybrids",
        subtitle = "All annotated genes",
        x = "Strain", y = "Gene")

ggsave("canme_all_mutations.pdf", mutated_gene_plot,
    width = 5, height = 14, limitsize = FALSE)

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

ggsave("canme_selected_mutations.pdf", select_gene_plot,
    width = 6, height = 7, limitsize = FALSE)

# Extract data for Vlad's lollipop plot
lp_data <- molten_data %>%
    # Select genes of interest
    filter(em_Preferred_name %in% c("FKS1", "MRR1")) %>%
    # Mark inherited/new/reverted SNPs in switch conditions
    group_by(X.Uploaded_variation) %>%
    mutate(switch = case_when(
        eStrain %in% switch_eStrains & any(eStrain %in% names(switch_eStrains)) ~ "inherited",
        eStrain %in% names(switch_eStrains) & all(!(eStrain %in% switch_eStrains)) ~ "reverted",
        eStrain %in% switch_eStrains ~ "new",
        eStrain %in% names(switch_eStrains) ~ "new",
        .default = "no_switch"
    )) %>%
    ungroup() %>%
    # Filter out inherited mutations
    # because we want to see when the snp originated
    filter(switch != "inherited") %>%
    # Remove those that did not originate in the expected regime
    # FKS1 is related to ANI resistance, so we don't expect it in F/AinF
    # MRR1 is related to FLZ resistance, so we don't expect it in A/FinA
    filter(
        !(em_Preferred_name == "FKS1" & cond %in% c("F", "AinF"))
    ) %>%
    filter(
        !(em_Preferred_name == "MRR1" & cond %in% c("A", "FinA"))
    ) %>%
    # Extract relevant columns
    select(em_Preferred_name, Protein_position, Amino_acids, Consequence) %>%
    # Cluster to unique snps in any strain
    unique()

# Save
write.csv(file = "canme_lollipop_data.csv", lp_data, row.names = FALSE)