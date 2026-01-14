#!/usr/bin/env Rscript
library(tidyverse)

# All LOH must be called against the same reference genome
args <- commandArgs(trailingOnly = TRUE)
reference_id <- "s172"

# Get files
loh_blocks <- as.list(list.files(pattern = "LOH_blocks.tsv"))
het_blocks <- as.list(list.files(pattern = "exp.het_blocks.bed"))
new_blocks <- list.files(pattern = "wo.*.bed")
new_blocks <- as.list(new_blocks[file.size(new_blocks) > 0])
genome_fai <- list.files(pattern = ".fai")

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
    filter(V3 - V2 > 2000)

# Read genome index
# Filter large scaffolds (for low-quality assemblies)
genome_idx <- read.table(genome_fai) %>%
    mutate(V1 = case_when(
        V1 == "NC_018292.1" ~ "1",
        V1 == "NC_018295.1" ~ "2",
        V1 == "NC_018296.1" ~ "3",
        V1 == "NC_018297.1" ~ "4",
        V1 == "NC_018298.1" ~ "5",
        V1 == "NC_018300.1" ~ "6",
        V1 == "NC_018301.1" ~ "7",
        V1 == "NC_018302.1" ~ "8",
        .default = V1
    ))

largest_scaffolds <- genome_idx %>%
    arrange(desc(V2)) %>%
    select(V1) %>%
    unique() %>%
    slice_head(n = 20) %>%
    pull(V1) %>%
    intersect(., plot_data$V1)

plot_data <- plot_data %>%
    mutate(sample_id = case_when(
        sample_id == "EPOM_037_Y10A" ~ "Co_YPD_10A",
        sample_id == "EPOM_087_A10E" ~ "Co_ANI_10E",
        sample_id == "EPOM_088_A10G" ~ "Co_ANI_10G",
        sample_id == "EPOM_133_F10A" ~ "Co_FLZ_10A",
        sample_id == "EPOM_134_F10C" ~ "Co_FLZ_10C",
        sample_id == "EPOM_183_AF10E" ~ "Co_ANIFLZ_10E",
        sample_id == "EPOM_229_FinA10A" ~ "Co_FinA_10A",
        sample_id == "EPOM_230_FinA10C" ~ "Co_FinA_10C",
        sample_id == "EPOM_280_AinF10G" ~ "Co_AinF_10G",
        sample_id == "EPOM_T0_10" ~ "Co_T0_10",
        .default = sample_id
    )) %>%
    mutate(sample_id = factor(sample_id,
        levels = c("Co_T0_10", "Co_YPD_10A", "Co_ANI_10E", "Co_ANI_10G",
            "Co_AinF_10G", "Co_FLZ_10A", "Co_FinA_10A", "Co_FLZ_10C", 
            "Co_FinA_10C", "Co_ANIFLZ_10E")
    )) %>%
    mutate(V1 = case_when(
        V1 == "NC_018292.1" ~ "1",
        V1 == "NC_018295.1" ~ "2",
        V1 == "NC_018296.1" ~ "3",
        V1 == "NC_018297.1" ~ "4",
        V1 == "NC_018298.1" ~ "5",
        V1 == "NC_018300.1" ~ "6",
        V1 == "NC_018301.1" ~ "7",
        V1 == "NC_018302.1" ~ "8",
        .default = V1
    ))
# Plot LOH across species for each scaffold/chromosome
plot <- ggplot(plot_data, aes(
        xmin = V2, xmax = V3, ymin = 0, ymax = 1
        )) +
    geom_rect(data = merge(genome_idx, as.data.frame(unique(plot_data$sample_id))),
        aes(xmin = 1, xmax = V2),
        fill = "white") +
    geom_rect(aes(fill = V10), colour = NA) +
    scale_fill_manual(name = "LOH",
        values = c("#5cb8e6ff", "#df3d3dff", "#ffe69e"),
        limits = c("REF", "ALT", "HET"),
        labels = c("LOH to A", "LOH to B", "Heterozygous (no LOH)")
    ) +
    scale_x_continuous(expand = c(0.05,0)) +
    facet_grid(sample_id ~ V1,
        scales = "free", space = "free_x") +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "bottom",
        panel.spacing = unit(0, "lines"))

ggsave(paste0(reference_id, "_jloh_plot.pdf"),
    plot,
    height = 5,
    width = 5,
    limitsize = FALSE)

