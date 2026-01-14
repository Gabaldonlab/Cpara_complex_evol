# Get sample names
source "scripts/utils/sample_names.sh"
source "scripts/utils/set_paths.sh"

module load bedtools

# Intersect the evolved CNV with the GFF for CANME
for (( i=0; i<${#ALL_CANME_STRAINS[@]}; i++ )); do
    strain="${ALL_CANME_STRAINS[$i]}"
    echo $strain
    # -wa retrieves the entire gene in the GFF if it has any overlap
    bedtools intersect -wa \
        -a data/references/CANME_BP57_mit_BP57.gff \
        -b $OUT_PATH/${strain}/final_CNVcalling.final.bed > \
        $ANNOT_PATH/${strain}.CNV_genes.gff
    # default retrieves only the overlapping region
    bedtools intersect \
        -a data/references/CANME_BP57_mit_BP57.gff \
        -b $OUT_PATH/${strain}/final_CNVcalling.final.bed > \
        $ANNOT_PATH/${strain}.CNV_regions.gff
    # get the overlapping region from the CNV tab
    bedtools intersect \
        -a $OUT_PATH/${strain}/final_CNVcalling.final.bed \
        -b $ANNOT_PATH/${strain}.CNV_regions.gff |
        sort -k1,1 -k2,2n > \
        $ANNOT_PATH/${strain}.CNV_regions.bed
    # -wo prints both entries as in a sort of merge
    bedtools intersect -wo \
        -a $OUT_PATH/${strain}/final_CNVcalling.final.bed  \
        -b data/references/CANME_BP57_mit_BP57.gff | \
        grep $'\tmRNA\t' > \
        $ANNOT_PATH/${strain}.CNV_GFF_intersect.bed
    # Add strain to combine later
    awk -v strain="${strain}" \
        'BEGIN {OFS="\t"} {print $0, strain}' \
        $ANNOT_PATH/${strain}.CNV_GFF_intersect.bed > \
        $ANNOT_PATH/${strain}.CNV_GFF_merge.bed
done

# Intersect the evolved CNV with the GFF for CANOR
for (( i=0; i<${#ALL_CANOR_STRAINS[@]}; i++ )); do
    strain="${ALL_CANOR_STRAINS[$i]}"
    echo $strain
    # -wa retrieves the entire gene in the GFF if it has any overlap
    bedtools intersect -wa \
        -a data/references/CANOR_CO90125_mit_MCO471.gff \
        -b $OUT_PATH/${strain}/final_CNVcalling.final.bed > \
        $ANNOT_PATH/${strain}.CNV_genes.gff
    # default retrieves only the overlapping region
    bedtools intersect \
        -a data/references/CANOR_CO90125_mit_MCO471.gff \
        -b $OUT_PATH/${strain}/final_CNVcalling.final.bed > \
        $ANNOT_PATH/${strain}.CNV_regions.gff
    # get the overlapping region from the CNV tab
    bedtools intersect \
        -a $OUT_PATH/${strain}/final_CNVcalling.final.bed \
        -b $ANNOT_PATH/${strain}.CNV_regions.gff |
        sort -k1,1 -k2,2n > \
        $ANNOT_PATH/${strain}.CNV_regions.bed
    # -wo prints both entries as in a sort of merge
    bedtools intersect -wo \
        -a $OUT_PATH/${strain}/final_CNVcalling.final.bed  \
        -b data/references/CANOR_CO90125_mit_MCO471.gff | \
        grep $'\tmRNA\t' > \
        $ANNOT_PATH/${strain}.CNV_GFF_intersect.bed
    # Add strain to combine later
    awk -v strain="${strain}" \
        'BEGIN {OFS="\t"} {print $0, strain}' \
        $ANNOT_PATH/${strain}.CNV_GFF_intersect.bed > \
        $ANNOT_PATH/${strain}.CNV_GFF_merge.bed
done