# Set environment
library(tidyverse)
library(data.table)
setwd("/gpfs/projects/bsc40/current/aredondo/evolved_hybrids/EPOM_workflow")
source("scripts/workflow/templates/canor_strains.r")

gff_file <- unlist(list.files(
    "data/references",
    recursive = TRUE,
    full.names = TRUE,
    pattern = "CANOR.*gff"
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

    # Select annotation tables and vcf files
    file_vec <- list.files(
        "results/annot_snps",
        recursive = TRUE,
        full.names = TRUE,
        pattern = "annotated")
    vcf_files <- list.files(
        "results/variant_subset",
        recursive = TRUE,
        full.names = TRUE,
        pattern = "wo_.*Y.*.vcf.gz$"        
    )
    vcfinloh_files <- list.files(
        "results/variant_loh",
        recursive = TRUE,
        full.names = TRUE,
        pattern = "wo_.*in_.*.vcf$"        
    )
    # Select files
    file_match <- unlist(sapply(
        c(clade_vec, t0_strains[clade], ypd_strains[clade]),
        FUN = function(name){grep(name, file_vec)}
    ))
    vcf_match <- unlist(sapply(
        c(clade_vec, t0_strains[clade], ypd_strains[clade]),
        FUN = function(name){grep(name, vcf_files)}
    ))
    vcfinloh_match <- unlist(sapply(
        c(clade_vec, t0_strains[clade], ypd_strains[clade]),
        FUN = function(name){grep(name, vcfinloh_files)}
    ))
    file_vec <- file_vec[file_match]
    vcf_files <- unique(vcf_files[vcf_match])
    vcfinloh_files <- unique(vcfinloh_files[vcfinloh_match])
    vcfinloh_files <- vcfinloh_files[file.size(vcfinloh_files) != 0]
    # Load annotation tables
    annot_list <- sapply(file_vec, FUN = function(file){
        read.table(file, sep = "\t", header = TRUE,
            comment.char = "", colClasses = c(Gene = "character"),
    )
    }, simplify = FALSE)
    names(annot_list) <- gsub("_annotated_variants.tab", "", gsub(".*/", "", names(annot_list)))
    # Load vcf files
    vcf_list <- sapply(vcf_files, FUN = function(file){
        fread(file, sep = "\t", header = TRUE) %>%
            rowwise() %>%
            mutate(GT = strsplit(SAMPLE, ":")[[1]][1])
    }, simplify = FALSE)
    names(vcf_list) <- gsub("_variants.*", "", gsub(".*/", "", names(vcf_list)))
    vcf_data <- vcf_list %>% bind_rows(.id = "eStrain")
    # Load vcfinloh files
    print(vcfinloh_files)
    if (length(vcfinloh_files) != 0) { 
        vcfinloh_list <- sapply(vcfinloh_files, FUN = function(file){
            fread(file, sep = "\t", header = FALSE, col.names =
                c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMPLE")) %>%
                rowwise() %>%
                mutate(GT = strsplit(SAMPLE, ":")[[1]][1])
        }, simplify = FALSE)
        names(vcfinloh_list) <- gsub("_variants.*", "", gsub(".*/", "", names(vcfinloh_list)))
        vcfinloh_data <- vcfinloh_list %>% bind_rows(.id = "eStrain") %>% mutate(inloh = "inloh")
    } else {
        vcfinloh_data <- data.frame(ID = character(), GT = character(), eStrain = character(), inloh = character())
    }
    # Prepare data for plotting
    molten_data <- bind_rows(annot_list, .id = "eStrain") %>%
        mutate(snp = paste0(Protein_position, Amino_acids)) %>%
        mutate(eStrain = factor(eStrain,
            levels = order_eStrain)) %>%
        mutate(cond = conditions_canor[order_eStrain[eStrain]]) %>%
        mutate(strain = clade) %>%
        mutate(tr_id = paste0("rna-", Feature)) %>%
        # Add genotype from the vcfs
        left_join(vcf_data %>% select(ID, GT, eStrain),
            by = join_by(X.Uploaded_variation == ID, eStrain == eStrain)) %>%
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
        )) %>%
        # Add tag to the SNPs in LOH
        left_join(vcfinloh_data %>% select(ID, eStrain, inloh),
            by = join_by(X.Uploaded_variation == ID, eStrain == eStrain))
    print(head(molten_data))
    # Save the data
    molten_snp_list[[clade]] <- molten_data
}

# Merge data from all clades
molten_snp_all <- bind_rows(molten_snp_list)
# Save
write.csv(
    molten_snp_all,
    "results/annot_plots/snp_data_canor.csv",
    row.names = FALSE
)