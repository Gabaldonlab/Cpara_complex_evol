// Process - PerSVade annotate_small_vars
process PSV_ANNOTSNPS {
    publishDir "${params.output_path}/annot_snps", mode: "copy"
    conda params.conda.psv_env

    executor "slurm"
    clusterOptions "-q ${params.queue} -c ${params.psv_annotsnps.threads}"
    time params.psv_annotsnps.time
    array params.psv_annotsnps.array_size
    cpus params.psv_annotsnps.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(reference_file), \
        path(index), \
        path(gff), \
        val(mitochondria), \
        path(vcf)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}_annotated_variants.tab")

    script:
    """
    echo "Starting PSV_ANNOTSNPS" `date +%H-%M-%S`
    python ${params.psv_exe} \
        annotate_small_vars \
        -o ${sample_id} \
        -r ${reference_file} \
        --mitochondrial_chromosome ${mitochondria} \
        --merged_vcf ${vcf} \
        -gff ${gff} \
        -mcode 3 -gcode 12 \
        --threads ${task.cpus} \
        --fraction_available_mem 0.95
   
    for file in \$(ls ${sample_id})
    do
    [ -f ${sample_id}/"\$file" ] && mv ${sample_id}/"\$file" ${sample_id}_\${file}
    done

    echo "Finished PSV_ANNOTSNPS" `date +%H-%M-%S`
    """
}

// Process - PerSVade annotate_SVs
process PSV_ANNOTSVS {
    publishDir "${params.output_path}/annot_svs", mode: "copy"
    conda params.conda.psv_env

    executor "slurm"
    clusterOptions "-q ${params.queue} -c ${params.psv_annotsnps.threads}"
    time params.psv_annotsnps.time
    array params.psv_annotsnps.array_size
    cpus params.psv_annotsnps.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(reference_file), \
        path(index), \
        path(gff), \
        val(mitochondria), \
        path(vcf, stageAs: "SV_and_CNV_variant_calling.vcf")

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}*.tab")

    script:
    """
    echo "Starting PSV_ANNOTSVS" `date +%H-%M-%S`
    python ${params.psv_exe} \
        annotate_SVs \
        -o ${sample_id} \
        -r ${reference_file} \
        --mitochondrial_chromosome ${mitochondria} \
        -gff ${gff} \
        --SV_CNV_vcf ${vcf} \
        -mcode 3 -gcode 12 \
        --threads ${task.cpus} \
        --fraction_available_mem 0.95

    for file in \$(ls ${sample_id})
    do
    [ -f ${sample_id}/"\$file" ] && mv ${sample_id}/"\$file" ${sample_id}_\${file}
    done

    echo "Finished PSV_ANNOTSVS" `date +%H-%M-%S`
    """
}

// Process - Merge GFF with PSV annotations
process MERGE_GFFSNP {
    publishDir "${params.output_path}/annot_plots", mode: "move"

    executor "local"

    input:
    tuple \
        val(ref_id), \
        path(annots), \
        path(gff)

    output:
    tuple \
        val(ref_id), \
        path("*.csv")

    script:
    """
    echo "Starting MERGE_GFFSNP" `date +%H-%M-%S`
    if [[ "${ref_id}" == "CANOR_CO90125_mit_MCO471" ]]
    then
        Rscript ${projectDir}/templates/snp_annot_canor.r \
            ${projectDir}/templates/canor_strains.r
    elif [[ "${ref_id}" == "CANME_BP57_mit_BP57" ]]
    then
        Rscript ${projectDir}/templates/snp_annot_canme.r \
            ${projectDir}/templates/canme_strains.r
    fi
    echo "Starting MERGE_GFFSNP" `date +%H-%M-%S`
    """

}

// Process - Plot SNPs
process SNP_PLOT {
    publishDir "${params.output_path}/annot_plots", mode: "move"

    executor "local"

    input:
    tuple \
        val(ref_id), \
        path(snp_csv)

    output:
    path("*.pdf")
    path("*.csv")

    script:
    """
    echo "Starting SNP_PLOT" `date +%H-%M-%S`
    if [[ "${ref_id}" == "CANOR_CO90125_mit_MCO471" ]]
    then
        Rscript ${projectDir}/templates/snp_plot_canor.r \
            ${projectDir}/templates/canor_strains.r
    elif [[ "${ref_id}" == "CANME_BP57_mit_BP57" ]]
    then
        Rscript ${projectDir}/templates/snp_plot_canme.r \
            ${projectDir}/templates/canme_strains.r
    fi
    echo "Finished SNP_PLOT" `date +%H-%M-%S`
    """
}

// Process - Merge GFF with LOH blocks
process MERGE_GFFLOH {
    publishDir "${params.output_path}/annot_loh", mode: "copy"

    executor "local"

    input:
    tuple \
        val(sample_id), \
        val(ref_id), \
        path(gff), \
        path(bed)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}.LOH_genes.gff"), \
        path("${sample_id}.LOH_region.gff"), \
        path("${sample_id}.LOH_regions.bed")

    script:
    """
    echo "Starting MERGE_GFFLOH" `date +%H-%M-%S`
    # -wa retrieves the entire gene in the GFF if it has any overlap
    bedtools intersect -wa \
        -a ${gff} \
        -b ${bed} > \
        ${sample_id}.LOH_genes.gff
    # default retrieves only the overlapping region
    bedtools intersect \
        -a ${gff} \
        -b ${bed} > \
        ${sample_id}.LOH_region.gff
    # get the overlapping blocks from the LOH tab with the allele info
    bedtools intersect \
        -a ${bed} \
        -b ${gff} > \
        ${sample_id}.LOH_regions.bed
    echo "Finished MERGE_GFFLOH" `date +%H-%M-%S`
    """
}