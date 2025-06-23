# Libraries
library(ggpubr)
library(tidyverse)
library(data.table)

# Parameters for the filtering
args <- commandArgs(trailingOnly = TRUE)
max_del_cov <- args[1]
min_dup_cov <- args[2]

# Load the CNV regions
# Contains the number of copies
file_vec <- list.files(
    ".", 
    recursive = TRUE,
    full.names = TRUE,
    pattern = "unfiltered_final_CNVcalling.tab")

cnv_tab_list <- sapply(file_vec, FUN = function(file){
    print(file)
    fread(file)
}, simplify = FALSE)

cnv_filt_list <- lapply(cnv_tab_list, FUN = function(tab){
    tab <- tab %>%
        mutate(cov_filt_pass = case_when(
            SVTYPE == "DEL" ~ median_coverage_corrected <= max_del_cov,
            SVTYPE == "DUP" ~ median_coverage_corrected >= min_dup_cov,
            .default = TRUE
        )) %>%
        filter(cov_filt_pass) %>%
        select(-cov_filt_pass)
})

write.table(
    cnv_filt_list[[1]], 
    paste0("final_CNVcalling.tab"),
    sep = "\t",
    col.names = TRUE,
    row.names = FALSE,
    quote = FALSE)