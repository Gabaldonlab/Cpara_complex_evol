# Set environment
library(tidyverse)

# Source the script to set the environment
args <- commandArgs(trailingOnly = TRUE)
source(args[1])

gff_file <- unlist(list.files(
    recursive = TRUE,
    full.names = TRUE,
    pattern = "gff"
))

# Load GFF
gff <- read.table(
    gff_file,
    header = FALSE, sep = "\t", as.is = c(rep(FALSE, 8), TRUE),
    col.names = c("Scaffold", "Source", "Feature_type", "Start", "End",
        "Score", "Strand", "Frame", "Misc")
)
# Retain only CDS features (the ones with annotations)
gff <- gff %>% filter(Feature_type == "mRNA")

# Parse GFF to extract the gene ID and em_Preferred_name from the Misc field
# Split the Misc field in its subfields
misc <- lapply(strsplit(gff$Misc, ";"), strsplit, split = "=")
misc <- lapply(misc, FUN = function(list){
    names(list) <- sapply(list, "[", 1)
    list
})
# Extract the ID from the Misc field
gff$Feature <- unlist(lapply(lapply(misc, "[[", "ID"), "[", 2))
# Extract the em_Preferred_name of Misc field
epn <- lapply(lapply(misc, "[[", "em_Preferred_name"), "[", 2)
# Remove nulls
print(epn)
epn[!(sapply(epn, is.character))] <- NA
epn <- unlist(epn)
gff$em_Preferred_name <- epn

# Add length variable
gff$Length <- gff$End - gff$Start

# Select annotation tables
file_vec <- list.files(
    recursive = TRUE,
    full.names = TRUE,
    pattern = "annotated"
)
# Load annotation tables
annot_list <- sapply(file_vec, FUN = function(file){
    cat("Reading", file, "\n")
    read.table(file, sep = "\t", header = TRUE, comment.char = "")
}, simplify = FALSE)

names(annot_list) <- sapply(names(annot_list), FUN = function(name){
    gsub("_annotated_variants.tab", "", gsub("\\./", "", name))
})

cat("Saving data\n")
# Save the SNP data in a tsv for easier reuse
snp_data_canme <- bind_rows(annot_list, .id = "eStrain") %>%
        mutate(snp = paste0(Protein_position, Amino_acids)) %>%
        mutate(cond = conditions_canme[eStrain]) %>%
        # Add Preferred_name from the GFF
        left_join(gff %>% select(Feature, em_Preferred_name),
            by = join_by(Feature == Feature)) %>%
        # Filter only protein-altering mutations
        filter(!(Consequence %in% c(
            "upstream_gene_variant",
            "downstream_gene_variant",
            "synonymous_variant"))) %>%
        # Replace Preferred_name with Feature when empty
        mutate(em_Preferred_name = case_when(
             em_Preferred_name == "-" | is.na(em_Preferred_name) ~ Feature,
             .default = em_Preferred_name
        ))

write.csv(
    snp_data_canme,
    "snp_data_canme.csv",
    row.names = FALSE
)
