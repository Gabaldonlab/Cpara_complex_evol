# Set environment
library(tidyverse)
setwd("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow")
source("scripts/workflow/templates/canor_strains.r")

gff_file <- "data/references/CANOR_CO90125_mit_MCO471.gff"

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
gff$Feature <- unlist(lapply(lapply(misc, "[[", "Name"), "[", 2))
# Extract the em_Preferred_name of Misc field
epn <- lapply(lapply(misc, "[[", "em_Preferred_name"), "[", 2)
# Remove nulls
epn[!(sapply(epn, is.character))] <- NA
epn <- unlist(epn)
gff$em_Preferred_name <- epn

# Add length variable
gff$Length <- gff$End - gff$Start

# Select annotation tables
file_vec <- list.files(
    "results/annot_svs",
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
    gsub("_annotated_variants.tab", "", gsub(".*/", "", name))
})

annot_list <- annot_list[order_eStrain_canor]

cat("Saving data\n")
# Save the SVs data in a tsv for easier reuse
svs_data_canor <- bind_rows(annot_list, .id = "eStrain") %>%
        mutate(cond = conditions_canor[eStrain]) %>%
        # Add Preferred_name from the GFF
        left_join(gff %>% select(Feature, em_Preferred_name),
            by = join_by(Feature == Feature)) %>%
        # Filter only protein-altering mutations
        filter(!(Consequence %in% c(
            "upstream_gene_variant",
            "downstream_gene_variant",
            "upstream_gene_variant_BND",
            "downstream_gene_variant_BND",
            "synonymous_variant"))) %>%
        # Replace Preferred_name with Feature when empty
        mutate(em_Preferred_name = case_when(
             em_Preferred_name == "-" | is.na(em_Preferred_name) ~ Feature,
             .default = em_Preferred_name
        ))

write.csv(
    svs_data_canor,
    "results/annot_plots/svs_data_canor.csv",
    row.names = FALSE
)
