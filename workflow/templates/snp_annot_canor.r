# Set environment
library(tidyverse)

# Source the script to set the environment
args <- commandArgs(trailingOnly = TRUE)
source(args[1])

gff_file <- unlist(list.files(
    ".",
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

# Parse GFF to extract the gene ID and locus_tag from the Misc field
# Split the Misc field in its subfields
misc <- lapply(strsplit(gff$Misc, ";"), strsplit, split = "=")
misc <- lapply(misc, FUN = function(list){
    names(list) <- sapply(list, "[", 1)
    list
})
# Extract the ID from the Misc field
gff$ID <- unlist(lapply(lapply(misc, "[[", "ID"), "[", 2))

# Extract the em_Preferred_name of Misc field
epn <- lapply(lapply(misc, "[[", "em_Preferred_name"), "[", 2)
# Remove nulls
epn[!(sapply(epn, is.character))] <- NA
epn <- unlist(epn)
gff$em_Preferred_name <- epn

# Extract the locus_tag of Misc field
lts <- lapply(lapply(misc, "[[", "locus_tag"), "[", 2)
# Remove nulls
lts[!(sapply(lts, is.character))] <- "locus_tag:empty"
lts <- unlist(lts)
gff$locus_tag <- lts
# Add length variable
gff$Length <- gff$End - gff$Start

# Set variables to iterate over clades
clades = c(
    "pa_strains" = "pa_strains",
    "c1_strains" = "c1_strains",
    "c2_strains" = "c2_strains",
    "c3_strains" = "c3_strains",
    "c4_strains" = "c4_strains",
    "pb_strains" = "pb_strains"
)

strain_names <- list(
    "pa_strains" = pa_strains,
    "c1_strains" = c1_strains,
    "c2_strains" = c2_strains,
    "c3_strains" = c3_strains,
    "c4_strains" = c4_strains,
    "pb_strains" = pb_strains
)

molten_snp_list <- list()
for (clade in clades){
    print(clade)
    # Choose clade
    clade_vec <- strain_names[[clade]]

    # Select annotation tables in the clade
    file_vec <- list.files(
        ".",
        recursive = TRUE,
        full.names = TRUE,
        pattern = "annotated")
    file_match <- unlist(sapply(
        c(clade_vec, t0_strains[clade], ypd_strains[clade]),
        FUN = function(name){grep(name, file_vec)}
    ))
    file_vec <- file_vec[file_match]

    # Load annotation tables
    annot_list <- sapply(file_vec, FUN = function(file){
        read.table(file, sep = "\t", header = TRUE,
            comment.char = "", colClasses = c(Gene = "character"),
    )
    }, simplify = FALSE)
    names(annot_list) <- gsub("_annotated_variants.tab", "", gsub("\\./", "", names(annot_list)))

    # Prepare data for plotting
    molten_data <- bind_rows(annot_list, .id = "eStrain") %>%
        mutate(snp = paste0(Protein_position, Amino_acids)) %>%
        mutate(eStrain = factor(eStrain,
            levels = order_eStrain)) %>%
        mutate(cond = conditions_canor[order_eStrain[eStrain]]) %>%
        mutate(strain = clade) %>%
        mutate(tr_id = paste0("rna-", Feature)) %>%
        # Add Preferred_name and locus_tag from the GFF
        left_join(gff %>% select(ID, em_Preferred_name),
            by = join_by(tr_id == ID)) %>%
        left_join(gff %>% select(ID, locus_tag),
            by = join_by(tr_id == ID)) %>%
        # Filter only protein-altering mutations
        filter(!(Consequence %in% c(
            "upstream_gene_variant",
            "downstream_gene_variant",
            "synonymous_variant"))) %>%
        # Replace Preferred_name with locus_tag when empty
        mutate(em_Preferred_name = case_when(
            em_Preferred_name == "-" | is.na(em_Preferred_name) ~ locus_tag,
            .default = em_Preferred_name
        ))

    # Save the data
    molten_snp_list[[clade]] <- molten_data
}

# Merge data from all clades
molten_snp_all <- bind_rows(molten_snp_list)
# Save
write.csv(
    molten_snp_all,
    "snp_data_canor.csv",
    row.names = FALSE
)