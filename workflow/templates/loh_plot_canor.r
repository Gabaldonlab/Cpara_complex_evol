# Set environment
library(tidyverse)
setwd("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow")
source("scripts/workflow/templates/canor_strains.r")

# Load data
molten_loh_all <- read.csv("results/annot_plots/loh_data_canor.csv")

molten_loh_all <- molten_loh_all %>%
    group_by(eStrain, geneID) %>%
    mutate(total_region_length = sum(region_length)) %>%
    ungroup() %>% 
    # Filter out the strains with whole genome LOH
    filter(!(eStrain %in% c("EPOM_088_A10G", "EPOM_133_F10A", "EPOM_183_AF10E",
        "EPOM_229_FinA10A", "EPOM_280_AinF10G"))) %>%
    # Filter out parentals
    filter(!(strain %in% c("Co90125", "SY36")))

# Plot
plot_loh_all <- ggplot(molten_loh_all,
    #filter(total_region_length / gene_length > 0.8),
    aes(
        x = eStrain, y = reorder(Preferred_name, rev(start)),
        alpha = total_region_length / gene_length,
        fill = allele
    )) +
    geom_tile(colour = "gray70") +
    facet_grid(chromosome ~ strain, scales = "free", space = "free") +
    theme_bw() +
    theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1),
        strip.background = element_rect(linetype = "blank"),
        strip.text.y = element_text(angle = 0),
        strip.clip = "off") +
    labs(title = "Overlap of genes with LOH blocks",
        subtitle = "Corrected for existing LOH at T0",
        caption = "LOH blocks in T0 strains were removed\nLOH blocks were not filtered by size\nGenes with any coverage of LOH were included",
        x = "Evolved strain",
        y = "Gene",
        alpha = "Extent covered by LOH")

ggsave("results/annot_plots/canor_loh_all_genes.pdf",
    plot_loh_all,
    width = 25,
    height = 200,
    limitsize = FALSE)

# Plot selected genes
loh_selected_data <- molten_loh_all %>% 
        filter(Preferred_name %in% interest_genes) #%>% 
        #filter(total_region_length / gene_length > 0.5)

plot_loh_selected <- ggplot(loh_selected_data,
    aes(
        x = eStrain, y = reorder(Preferred_name, rev(start)),
        alpha = total_region_length / gene_length,
        fill = allele
    )) +
    geom_tile(colour = "gray70") +
    facet_grid(. ~ cond, scales = "free", space = "free") +
    theme_bw() +
    theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 60, vjust = 1, hjust = 1),
        strip.background = element_rect(linetype = "blank"),
        strip.text.y = element_text(angle = 0),
        strip.clip = "off") +
    labs(title = "Overlap of genes with LOH blocks",
        subtitle = "Corrected for existing LOH at T0",
        caption = "LOH blocks in T0 strains were removed\nLOH blocks were not filtered by size\nGenes with any coverage of LOH were included",
        x = "Evolved strain",
        y = "Gene",
        alpha = "Extent covered by LOH") +
    scale_y_discrete(limits = interest_genes)

ggsave("results/annot_plots/canor_loh_selected_genes.pdf",
    plot_loh_selected,
    width = 5,
    height = 8,
    limitsize = FALSE)

# Plot for the paper
loh_paperplot_data <- molten_loh_all %>% 
        filter(Preferred_name %in% interest_genes) %>% 
        #filter(total_region_length / gene_length > 0.5) %>%
        mutate(eStrain = gsub("EPOM_", "", eStrain))

plot_loh_paperplot <- ggplot(loh_paperplot_data,
    aes(
        x = eStrain,
        y = reorder(Preferred_name, rev(start)),
        fill = allele
    )) +
    geom_tile(colour = "gray70") +
    facet_grid(. ~ cond, scales = "free", space = "free") +
    theme_bw() +
    theme(legend.position = "bottom",
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
        strip.background = element_rect(linetype = "blank"),
        strip.text.y = element_text(angle = 0),
        strip.clip = "off") +
    labs(title = "Overlap of genes with LOH blocks",
        caption = "LOH blocks in T0 strains were removed\nLOH blocks were not filtered by size\nGenes with any coverage of LOH were included",
        x = "Evolved strain",
        y = "Gene")

ggsave("results/annot_plots/canor_loh_paperplot_genes.pdf",
    plot_loh_paperplot,
    width = 8,
    height = 5.2,
    limitsize = FALSE)