#!/usr/bin/env Rscript
library(tidyverse)
library(ggpubr)

# All LOH must be called against the same reference genome
args <- commandArgs(trailingOnly = TRUE)
reference_id <- "all_LOH"

# Source sample names
source("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow/scripts/workflow/templates/canor_strains.r")
source("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow/scripts/workflow/templates/canme_strains.r")

# Get files
loh_blocks <- as.list(list.files(pattern = "LOH_blocks.tsv"))
het_blocks <- as.list(list.files(pattern = "exp.het_blocks.bed"))
new_blocks <- list.files(pattern = "wo.*.bed")
new_blocks <- as.list(new_blocks[file.size(new_blocks) > 0])
canor_fai <- list.files("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow/data/references",
    pattern = "CANOR_CO90125_mit_MCO471.fna.fai",
    full.names = TRUE)
canme_fai <- list.files("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow/data/references",
    pattern = "CANME_BP57_mit_BP57.fna.fai",
    full.names = TRUE)

names(loh_blocks) <- sapply(loh_blocks,
    gsub, pattern = ".LOH_blocks.tsv", replacement = "")
names(het_blocks) <- sapply(het_blocks,
    gsub, pattern = ".exp.het_blocks.bed", replacement = "")
names(new_blocks) <- sapply(new_blocks,
    gsub, pattern = "_wo.*.bed", replacement = "")

# Read LOH data
loh_data <- lapply(loh_blocks, FUN = function(file){
    read.table(file, sep = "\t") %>%
        select(V1, V2, V3, V8, V10)
}) %>% bind_rows(.id = "sample_id")

het_data <- lapply(het_blocks, FUN = function(file){
    read.table(file, sep = "\t") %>%
        mutate(V10 = "HET") %>%
        mutate(V8 = "homo") %>%
        select(V1, V2, V3, V8, V10)
}) %>% bind_rows(.id = "sample_id")

new_data <- lapply(new_blocks, FUN = function(file){
    read.table(file, sep = "\t") %>%
        mutate(V10 = "NEW") %>%
        mutate(V8 = "homo") %>%
        select(V1, V2, V3, V8, V10)
}) %>% bind_rows(.id = "sample_id")

plot_data <- rbind(loh_data, het_data) %>%
    filter(V3 - V2 > 2000) %>%
    mutate(sample_fct = factor(sample_id,
        levels = c(canor_hybird_strains, order_allstrains_canme))
    )

# Define interest genes
interest_genes_canor <- data.frame(
    V1 = c("Co_Chr2", "Co_Chr3"),
    V2 = c(1597409, 399558),
    V3 = c(1603117, 402836)
)

interest_genes_canme <- data.frame(
    V1 = c("Cm_Chr1", "Cm_Chr1"),
    V2 = c(1960966, 2969495),
    V3 = c(1966701, 2969591)
)
# Read genome index
# Filter large scaffolds (for low-quality assemblies)
canor_idx <- read.table(canor_fai) %>%
    mutate(V1 = case_when(
        V1 == "NC_018292.1" ~ "Co_Chr1",
        V1 == "NC_018295.1" ~ "Co_Chr2",
        V1 == "NC_018296.1" ~ "Co_Chr3",
        V1 == "NC_018297.1" ~ "Co_Chr4",
        V1 == "NC_018298.1" ~ "Co_Chr5",
        V1 == "NC_018300.1" ~ "Co_Chr6",
        V1 == "NC_018301.1" ~ "Co_Chr7",
        V1 == "NC_018302.1" ~ "Co_Chr8",
        .default = V1
    ))

canme_idx <- read.table(canme_fai) %>%
    mutate(V1 = factor(case_when(
        V1 == "scaffold1" ~ "Cm_Chr1",
        V1 == "scaffold2" ~ "Cm_Chr2",
        V1 == "scaffold3" ~ "Cm_Chr3",
        V1 == "scaffold4" ~ "Cm_Chr4",
        V1 == "scaffold5" ~ "Cm_Chr5",
        V1 == "scaffold6" ~ "Cm_Chr6",
        V1 == "scaffold7" ~ "Cm_Chr7",
        V1 == "scaffold8" ~ "8",
        V1 == "scaffold9" ~ "9",
        .default = V1
    ), levels = c("Cm_Chr1", "Cm_Chr2", "Cm_Chr3", "Cm_Chr4", 
        "Cm_Chr5", "Cm_Chr6", "Cm_Chr7", "8", "9")))

plot_data <- plot_data %>%
    mutate(V1 = factor(case_when(
        V1 == "NC_018292.1" ~ "Co_Chr1",
        V1 == "NC_018295.1" ~ "Co_Chr2",
        V1 == "NC_018296.1" ~ "Co_Chr3",
        V1 == "NC_018297.1" ~ "Co_Chr4",
        V1 == "NC_018298.1" ~ "Co_Chr5",
        V1 == "NC_018300.1" ~ "Co_Chr6",
        V1 == "NC_018301.1" ~ "Co_Chr7",
        V1 == "NC_018302.1" ~ "Co_Chr8",
        V1 == "scaffold1" ~ "Cm_Chr1",
        V1 == "scaffold2" ~ "Cm_Chr2",
        V1 == "scaffold3" ~ "Cm_Chr3",
        V1 == "scaffold4" ~ "Cm_Chr4",
        V1 == "scaffold5" ~ "Cm_Chr5",
        V1 == "scaffold6" ~ "Cm_Chr6",
        V1 == "scaffold7" ~ "Cm_Chr7",
        V1 == "scaffold8" ~ "8",
        V1 == "scaffold9" ~ "9",
        .default = V1
    ), levels = c("Co_Chr1", "Co_Chr2", "Co_Chr3", "Co_Chr4", 
        "Co_Chr5", "Co_Chr6", "Co_Chr7", "Co_Chr8",
        "Cm_Chr1", "Cm_Chr2", "Cm_Chr3", "Cm_Chr4", 
        "Cm_Chr5", "Cm_Chr6", "Cm_Chr7", "8", "9"))) %>%
    mutate(species = case_when(
        sample_id %in% order_allstrains_canme ~ "CANME",
        sample_id %in% canor_hybird_strains ~ "CANOR",
        .default = NA
    ))
# Plot LOH across species for each scaffold/chromosome
plot_canor <- ggplot(plot_data %>% filter(species == "CANOR"), aes(
        xmin = V2, xmax = V3, ymin = 0, ymax = 1
        )) +
    geom_rect(data = merge(canor_idx, as.data.frame(unique(plot_data$sample_id))),
        aes(xmin = 1, xmax = V2),
        fill = "white") +
    geom_rect(aes(fill = V10), colour = NA) +
    scale_fill_manual(name = "LOH",
        values = c("#5cb8e6ff", "#df3d3dff", "#ffe69e"),
        limits = c("REF", "ALT", "HET"),
        labels = c("LOH to parental A", "LOH to parental B", "Heterozygous (no LOH)")
    ) +
    scale_x_continuous(expand = c(0.05,0)) +
    facet_grid(sample_fct ~ V1,
        scales = "free", space = "free_x") +
    theme_classic() +
    theme(plot.margin = margin(1,0.1,0.1,0.1, "cm"),
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0, "lines"))

plot_canme <- ggplot(plot_data %>% filter(species == "CANME"), aes(
        xmin = V2, xmax = V3, ymin = 0, ymax = 1
        )) +
    geom_rect(data = merge(canme_idx, as.data.frame(unique(plot_data$sample_id))),
        aes(xmin = 1, xmax = V2),
        fill = "white") +
    geom_rect(aes(fill = V10), colour = NA) +
    scale_fill_manual(name = "LOH",
        values = c("#4d804d", "#7f5f9c", "#ffe69e"),
        limits = c("REF", "ALT", "HET"),
        labels = c("LOH to Reference", "LOH to Alternative", "Heterozygous (no LOH)")
    ) +
    scale_x_continuous(expand = c(0.05,0)) +
    facet_grid(sample_fct ~ V1,
        scales = "free", space = "free_x") +
    theme_classic() +
    theme(plot.margin = margin(1,0.1,0.1,0.1, "cm"),
        strip.text.y = element_text(angle = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0, "lines"))


plot_combined <- ggarrange(plot_canor, plot_canme, ncol = 1,
    labels = c("A", "B"),
    heights = c(4 + 36, 4 + 16))

ggsave(paste0("../jloh_plot/", reference_id, "_jloh_plot.pdf"),
    plot_combined,
    height = 15,
    width = 15,
    limitsize = FALSE)

ggsave(paste0("../jloh_plot/", reference_id, "_jloh_plot.png"),
    plot_combined,
    height = 15,
    width = 15,
    limitsize = FALSE)

