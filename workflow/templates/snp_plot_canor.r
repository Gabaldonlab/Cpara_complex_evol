# Set environment
library(tidyverse)
setwd("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow")
source("scripts/workflow/templates/canor_strains.r")

# Load data
molten_snp_all <- read.csv("results/annot_plots/snp_data_canor.csv")

# Prepare data for plotting
plot_snp_data <- molten_snp_all %>%
        # Reorder factors
        mutate(cond = factor(cond, levels = c("ANI", "AinF", "ANIFLZ", "FLZ", "FinA"))) %>%
        # # Filter out SNPs in all eStrains of the same strain
        # group_by(strain) %>%
        # mutate(n_eStrain = length(unique(eStrain))) %>%
        # group_by(strain, snp) %>%
        # mutate(prevalence = n_distinct(eStrain)) %>%
        # filter(prevalence < n_eStrain) %>%
        # Mark inherited/new/reverted SNPs in switch conditions
        ungroup() %>%
        mutate(eStrain = factor(eStrain,
            levels = order_eStrain_canor)) %>%
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
        group_by(eStrain, tr_id) %>%
        mutate(cumm_snps = length(unique(X.Uploaded_variation))) %>%
        group_by(tr_id) %>%
        # Sort by number of strains with a mutation
        mutate(total_snps = length(unique(eStrain))) %>%
        arrange(total_snps) %>%
        # Modify the inloh categories
        mutate(inloh = case_when(
            eStrain %in% c("EPOM_133_F10A", "EPOM_144_F12G", "EPOM_183_AF10E", "EPOM_229_FinA10A") ~ "haploid",
            eStrain %in% c("EPOM_280_AinF10G", "EPOM_088_A10G") ~ "wgloh",
            grepl("12[A-Z]$", eStrain) ~ "parental",
            grepl("7[A-Z]$", eStrain) ~ "parental",
            .default = inloh
            )
        )

# Generate plot
full_mutated_gene_plot <- ggplot(
    plot_snp_data %>%
        mutate(em_Preferred_name = factor(em_Preferred_name,
            levels = unique(plot_snp_data$em_Preferred_name))),
    aes(x = eStrain,
        y = snp,
        fill = Consequence,
        label = cumm_snps,
        alpha = switch)) +
    geom_tile() +
    scale_alpha_manual(
            values = c("new" = 1, "inherited" = 0.1, "reverted" = 0.5),
            na.translate = FALSE
    ) +
    facet_grid(em_Preferred_name ~ cond, scales = "free", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom") +
    labs(title = "Missense variants in evolved hybrids",
        x = "Strain", y = "Gene")

ggsave("results/annot_plots/canor_all_mutations.pdf", full_mutated_gene_plot,
    width = 10, height = 17, limitsize = FALSE)

# Generate plot only with the selected genes
select_gene_plot <- ggplot(
    plot_snp_data %>% 
        # Filter out unnamed genes
        filter(em_Preferred_name != "") %>%
        # Subset selected genes
        filter(toupper(em_Preferred_name) %in% selected_genes),
    aes(x = eStrain,
        y = snp,
        fill = Consequence,
        colour = GT,
        alpha = switch)) +
    geom_tile() +
    geom_point(aes(shape = inloh), alpha = 1) +
    scale_alpha_manual(
        values = c("new" = 1, "inherited" = 0.2, "reverted" = 0.5),
        na.translate = FALSE
    ) +
    scale_colour_manual(
        values = c("0/1" = "white", "1/1" = "black"),
        na.translate = FALSE
    ) + 
    scale_shape_manual(
        values = c("parental" = "P", "haploid" = "H",
            "inloh" = "L", "wgloh" = "W"),
        na.translate = FALSE
    ) +
    scale_fill_manual(
        values = c("missense_variant" = "coral", "stop_gained" = "firebrick",
            "inframe_deletion" = "slateblue", "inframe_insertion" = "cyan",
            "frameshift_variant" = "goldenrod"),
        na.translate = FALSE
    ) +
    facet_grid(em_Preferred_name ~ cond, scales = "free", space = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        strip.text.y = element_text(angle = 0),
        legend.position = "bottom") +
    guides(fill = guide_legend(ncol = 1), colour = guide_legend(ncol = 1),
        shape = guide_legend(ncol = 1), alpha = guide_legend(ncol = 1) ) +
    labs(title = "Protein-altering variants in evolved hybrids",
        subtitle = "Selected genes",
        x = "Strain", y = "Mutation",
        fill = "Consequence", alpha = "Switch", shape = "Recombination",
        colour = "Genotype")

ggsave("results/annot_plots/canor_selected_mutations.pdf", select_gene_plot,
    width = 8, height = 10, limitsize = FALSE)

# Extract data for Vlad's lollipop plot
lp_data <- molten_snp_all %>%
    # Select genes of interest
    filter(em_Preferred_name %in% c("FKS1", "MRR1")) %>%
    # Mark inherited/new/reverted SNPs in switch conditions
    group_by(X.Uploaded_variation, eStrain) %>%
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
    select(eStrain, em_Preferred_name, Protein_position, Amino_acids, Consequence, GT) %>%
    # Cluster to unique snps in any strain
    unique()

# Save
write.csv(file = "results/annot_plots/canor_lollipop_data.csv", lp_data, row.names = FALSE)