// Process - JLOH stats
process JLOH_STATS {
    publishDir "${params.output_path}/jloh_stats"
    conda params.conda.jloh_env
   
    executor "slurm"
    clusterOptions "-q ${params.queue} -c ${params.jloh_stats.threads}"
    time params.jloh_stats.time
    array params.jloh_stats.array_size
    cpus params.jloh_stats.threads
    errorStrategy "ignore"
    
    input:
    tuple \
        val(sample_id), \
        path(variants)
    
    output:
    tuple \
        val(sample_id), \
        path("${sample_id}_jlohstats.txt")

    script:
    """
    echo "Starting JLOH_STATS" `date +%H-%M-%S`
    python ${params.jloh_exe} stats \
        --window-size ${params.jloh_stats.window_size} \
        --step-size ${params.jloh_stats.step_size} \
        --vcf ${variants} \
        --threads ${task.cpus} &> ${sample_id}_jlohstats.txt;
    echo "Finished JLOH_STATS" `date +%H-%M-%S`
    """
}

// Process - Parse JLOH stats
process PARSE_JLOHSTATS{
    executor "local"
    
    input:
    tuple \
        val(sample_id), \
        path(jlohstats)

    output:
    tuple \
        val(sample_id), \
        env(hetsnp), \
        env(homsnp)

    script:
    """
    echo "Starting PARSE_JLOHSTATS" `date +%H-%M-%S`
    hetsnp=\$(grep "50%" ${jlohstats} | cut -f 3)
    homsnp=\$(grep "50%" ${jlohstats} | cut -f 4)
    # Check if values are "nan" and replace them with 0
    [ "\$hetsnp" = "nan" ] && hetsnp=0
    [ "\$homsnp" = "nan" ] && homsnp=0
    echo "Finished PARSE_JLOHSTATS" `date +%H-%M-%S`
    """
}

// Process - JLOH extract
process JLOH_EXTRACT{
    publishDir "${params.output_path}/jloh_extract", mode: "copy"
    conda params.conda.jloh_env

    executor "slurm"
    clusterOptions "-q ${params.queue} -c ${params.jloh_extract.threads}"
    time params.jloh_extract.time
    array params.jloh_extract.array_size
    cpus params.jloh_extract.threads
    errorStrategy "ignore"

    input:
    // Tendrás que crear este canal mezclando cosas
    tuple \
        val(sample_id), \
        path(reference_file), \
        path(index), \
        path(alignment), \
        path(alignment_index), \
        path(variants), \
        val(hetsnp), \
        val(homsnp)

    output:
    tuple \
        val(sample_id), \
        path("*.LOH_blocks.bed"), \
        path("*.exp.het_blocks.bed"), \
        path("*.LOH_blocks.tsv"), \
        path("*.uncalled_blocks.bed")

    script:
    """
    echo "Starting JLOH_EXTRACT" `date +%H-%M-%S`
    python ${params.jloh_exe} extract \
        --vcf ${variants} \
        --bam ${alignment} \
        --ref ${reference_file} \
        --threads ${task.cpus} \
        --min-length ${params.jloh_extract.min_length} \
        --min-frac-cov ${params.jloh_extract.min_frac_cov} \
        --output-dir . \
        --sample ${sample_id} \
        --min-snps-kbp ${hetsnp},${homsnp}

    # Parse BED files to include the allele field
    awk -F "\\t" \
        'BEGIN {OFS="\\t"} NR > 1 { \$2=\$2-1; print \$1, \$2, \$3, \$10 }' \
        ${sample_id}.LOH_blocks.tsv > ${sample_id}.LOH_blocks.bed

    # Extract uncalled regions (not LOH or het)
    bedtools subtract \
        -a ${reference_file.simpleName}.bed \
        -b ${sample_id}.LOH_blocks.bed > \
        tmp.uncalled_blocks.bed

    bedtools subtract \
        -a tmp.uncalled_blocks.bed \
        -b ${sample_id}.exp.het_blocks.bed > \
        ${sample_id}.uncalled_blocks.bed

    rm tmp.uncalled_blocks.bed
    echo "Finished JLOH_EXTRACT" `date +%H-%M-%S`
    """
}

// Process - BEDtools LOH stats
process BED_LOH_STATS{
    publishDir "${params.output_path}/bed_loh_stats", mode: "copy"
    executor "local"

    input:
    tuple \
        val(sample_id), \
        path(loh_bed), \
        path(loh_tsv), \
        path(het_bed), \
        path(genome_file)

    output:
    tuple \
        val(sample_id), \
        path("*_stats_het_blocks.tsv"), \
        path("*_stats_LOH_blocks.tsv")

    script:
    """
    echo "Starting BED_LOH_STATS" `date +%H-%M-%S`
    bedtools summary -i ${het_bed} -g ${genome_file} > ${sample_id}_stats_het_blocks.tsv
    bedtools summary -i ${loh_bed} -g ${genome_file} > ${sample_id}_stats_LOH_blocks.tsv 
    echo "Finished BED_LOH_STATS" `date +%H-%M-%S`
    """
}

// Process - Merge BED files
process MERGE_BEDS{
    publishDir "${params.output_path}/jloh_extract", mode: "copy"
    executor "local"

    input:
    tuple \
        val(sample_id), \
        path(loh_bed), \
        path(het_bed), \
        path(loh_tsv), \
        path(uncalled_bed)

    output:
    tuple \
        val(sample_id), \
        path("*.LOH_blocks_and_uncalled.bed"), \
        path("*.exp.het_blocks.bed"), \
        path("*.LOH_blocks.tsv"), \
        path("*.uncalled_blocks.bed")

    script:
    """
    echo "Starting BED_INTERSECT" `date +%H-%M-%S`
    bedtools merge \
        -a ${loh_bed} -b ${uncalled_bed} > ${sample_id}.LOH_blocks_and_uncalled.bed
    echo "Finished BED_INTERSECT" `date +%H-%M-%S`
    """  
}

// Process - BEDtools intersect
process BED_INTERSECT{
    publishDir "${params.output_path}/jloh_extract", mode: "copy"
    executor "local"

    input:
    tuple \
        val(sample_id), \
        path(bed), \
        val(background_id), \
        path(bed_background)

    output:
    tuple \
        val(sample_id), \
        path("${bed.simpleName}_wo_${background_id}.bed")

    script:
    """
    echo "Starting BED_INTERSECT" `date +%H-%M-%S`
    bedtools subtract \
        -a ${bed} -b ${bed_background} > ${bed.simpleName}_wo_${background_id}.bed
    echo "Finished BED_INTERSECT" `date +%H-%M-%S`
    """    
}

// Process - JLOH plot
process JLOH_PLOT{
    publishDir "${params.output_path}/jloh_plot", mode: "copy"
    executor "local"

    input:
    tuple \
        val(reference_id), \
        path(all_jloh_files), \
        path(genome_index)

    output:
    path("${reference_id}_jloh_plot.pdf")

    script:
    """
    echo "Starting JLOH_PLOT" `date +%H-%M-%S`
    Rscript ${projectDir}/templates/jloh_plot_multiple.r ${reference_id}
    echo "Finished JLOH_PLOT" `date +%H-%M-%S`
    """     
}

// Process - JLOH plot
process JLOH_PLOT_VARS{
    publishDir "${params.output_path}/jloh_plot", mode: "copy"
    executor "local"

    input:
    tuple \
        val(reference_id), \
        path(all_files), \
        path(genome_index)

    output:
    path("${reference_id}_jloh_plot_w_vars.png")

    script:
    """
    echo "Starting JLOH_PLOT_VARS" `date +%H-%M-%S`
    Rscript ${projectDir}/templates/jloh_plot_w_vars.r ${reference_id}
    echo "Finished JLOH_PLOT_VARS" `date +%H-%M-%S`
    """     
}

// Process - BEDtools intersect with VCF
process BED_VCF_INTERSECT{
    publishDir "${params.output_path}/variant_loh", mode: "copy"
    executor "local"

    input:
    tuple \
        val(sample_id), \
        path(vcf), \
        path(bed)

    output:
    tuple \
        val(sample_id), \
        path("${vcf.simpleName}_in_${bed.simpleName}.vcf")

    script:
    """
    echo "Starting BED_VCF_INTERSECT" `date +%H-%M-%S`
    bedtools intersect \
        -a ${vcf} -b ${bed} > ${vcf.simpleName}_in_${bed.simpleName}.vcf
    echo "Finished BED_VCF_INTERSECT" `date +%H-%M-%S`
    """    
}