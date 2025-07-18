#!/usr/bin/env Rscript
library(tidyverse)

# All LOH must be called against the same reference genome
args <- commandArgs(trailingOnly = TRUE)
reference_id <- args[1]

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

plot_data <- rbind(loh_data, het_data)

# Read genome index
# Filter large scaffolds (for low-quality assemblies)
genome_idx <- read.table(genome_fai)

largest_scaffolds <- genome_idx %>%
    arrange(desc(V2)) %>%
    select(V1) %>%
    unique() %>%
    slice_head(n = 20) %>%
    pull(V1) %>%
    intersect(., plot_data$V1)

# Plot LOH across species for each scaffold/chromosome
plot <- ggplot(plot_data, aes(
        xmin = V2, xmax = V3, ymin = 0, ymax = 1
        )) +
    geom_rect(data = merge(genome_idx, as.data.frame(unique(plot_data$sample_id))),
        aes(xmin = 1, xmax = V2),
        fill = "white") +
    geom_rect(aes(fill = V10, alpha = V8)) +
    geom_rect(data = new_data, aes(ymin = 1, ymax = 1.5), fill = "slateblue") +
    scale_fill_manual(
        values = c("dodgerblue2", "firebrick", "#fffd87"),
        limits = c("REF", "ALT", "HET")
    ) +
    scale_alpha_manual(
        values = c(0.5, 1),
        limits = c("hemi", "homo")
    ) +
    facet_grid(sample_id ~ V1,
        scales = "free", space = "free_x") +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 0),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
    labs(title = "LOH blocks")

ggsave(paste0(reference_id, "_jloh_plot.pdf"),
    plot,
    height = 0.7 * min(20, length(unique(plot_data$sample_id))),
    width = 40,
    limitsize = FALSE)

