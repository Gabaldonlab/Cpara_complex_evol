#!/usr/bin/env Rscript
library(tidyverse)

# All LOH must be called against the same reference genome
args <- commandArgs(trailingOnly = TRUE)
reference_id <- args[1]

# Get files
loh_blocks <- as.list(list.files(pattern = ".LOH_blocks.tsv"))
het_blocks <- as.list(list.files(pattern = ".exp.het_blocks.bed"))
loh_alt <- as.list(list.files(pattern = "_lohAltVars.vcf"))
loh_ref <- as.list(list.files(pattern = "_lohRefVars.vcf"))
non_loh <- as.list(list.files(pattern = "_nonLohVars.vcf"))
new_hom <- as.list(list.files(pattern = "wo"))
genome_fai <- list.files(pattern = ".fai")

names(loh_blocks) <- sapply(loh_blocks,
    gsub, pattern = ".LOH_blocks.tsv", replacement = "")
names(het_blocks) <- sapply(het_blocks,
    gsub, pattern = ".exp.het_blocks.bed", replacement = "")
names(loh_alt) <- sapply(loh_alt,
    FUN = function(name){
        paste(strsplit(name, "_")[[1]][c(1,2,3)], collapse = "_")
    })
names(loh_ref) <- sapply(loh_ref,
    FUN = function(name){
        paste(strsplit(name, "_")[[1]][c(1,2,3)], collapse = "_")
    })
names(non_loh) <- sapply(non_loh,
    FUN = function(name){
        paste(strsplit(name, "_")[[1]][c(1,2,3)], collapse = "_")
    })
names(new_hom) <- sapply(new_hom,
    FUN = function(name){
        paste(strsplit(name, "_")[[1]][c(1,2,3)], collapse = "_")
    })

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

plot_data <- rbind(loh_data, het_data)

read_vcf <- function(file, class_name) {
    tryCatch(
        {
            print(file)
            read.table(file, sep = "\t") %>%
            select(V1, V2) %>%
            mutate(class = class_name)
        },
        error = function(cond) {
            print(paste("file", file, "seems empty"))
            data.frame(V1 = character(), V2 = integer())
        }
    )
}

loh_alt_data <- lapply(loh_alt, read_vcf, "het2hom") %>%
    bind_rows(.id = "sample_id")

loh_ref_data <- lapply(loh_ref, read_vcf, "het_lost") %>%
    bind_rows(.id = "sample_id")

non_loh_data <- lapply(non_loh, read_vcf, "het_snps") %>%
    bind_rows(.id = "sample_id")

new_hom_data <- lapply(new_hom, read_vcf, "new_hom") %>%
    bind_rows(.id = "sample_id")

var_data <- rbind(loh_alt_data, loh_ref_data, non_loh_data, new_hom_data) %>%
    mutate(class_y = case_when(
        class == "new_hom" ~ 4.5,
        class == "het2hom" ~ 4.0,
        class == "het_lost" ~ 3.5,
        class == "het_snps" ~ 3.0
    ))
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
plot <- ggplot() +
    geom_rect(data = merge(genome_idx, as.data.frame(unique(plot_data$sample_id))),
        aes(xmin = 1, xmax = V2, ymin = 1, ymax = 2),
        fill = "white") +
    geom_rect(data = plot_data,
        aes(xmin = V2, xmax = V3, ymin = 5, ymax = 6, fill = V10, alpha = V8)) +
    geom_point(data = var_data,
        aes(x = V2, y = class_y, colour = class),
        shape = 19, alpha = 0.3) +
    geom_density(data = var_data,
        aes(x = V2, y = after_stat(scaled), fill = class, colour = class),
        adjust = 0.25, alpha = 0.3, position = "stack") +
    scale_fill_manual(
        values = c("firebrick", "dodgerblue2", "#fffd87",
            "firebrick", "dodgerblue2", "#fffd87"),
        limits = c("REF", "ALT", "HET",
            "het_lost", "het2hom", "het_snps")
    ) +
    scale_colour_manual(
        values = c("dodgerblue4", "firebrick", "dodgerblue2", "#fffd87"),
        limits = c("new_hom", "het_lost", "het2hom", "het_snps")
    ) +
    scale_alpha_manual(
        values = c(0.5, 1),
        limits = c("hemi", "homo")
    ) +
    facet_grid(sample_id ~ V1,
        scales = "free", space = "free_x") +
    theme_classic() +
    theme(strip.text.y = element_text(angle = 0),
        #axis.text.y = element_blank(),
        #axis.ticks.y = element_blank(),
        axis.line.y = element_blank()) +
    labs(title = "LOH blocks with LOH variants")

ggsave(paste0(reference_id, "_jloh_plot_w_vars.png"),
    plot,
    height = 1 * min(20, length(unique(plot_data$sample_id))),
    width = 40,
    limitsize = FALSE)

