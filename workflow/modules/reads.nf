// Process - Filter raw reads
process PSV_TRIMQC {
    publishDir "${params.output_path}/reads_filtered_and_qc", mode: "copy"
    conda params.conda.psv_env

    executor "slurm"
    clusterOptions "-q gp_bscls -c ${params.psv_trimqc.threads}"
    time { params.psv_trimqc.time * task.attempt }
    array params.psv_trimqc.array_size
    cpus params.psv_trimqc.threads
    errorStrategy "retry"

    input:
    tuple \
        val(sample_id), \
        path(reads_fwd), \
        path(reads_rev)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}_trimmed_reads1.fastq.gz"), \
        path("${sample_id}_trimmed_reads2.fastq.gz"), \
        path("${sample_id}_trimmed_reads1_fastqc.zip"), \
        path("${sample_id}_trimmed_reads2_fastqc.zip"), \
        path("${sample_id}_trimmed_reads1_fastqcReport.html"), \
        path("${sample_id}_trimmed_reads2_fastqcReport.html")
    
    script:
    """
    echo "Starting PSV_TRIMQC" `date +%H-%M-%S`
    ${params.psv_exe} trim_reads_and_QC \
        -o ${sample_id} \
        -f1 ${reads_fwd} \
        -f2 ${reads_rev} \
        -thr ${task.cpus} \
        --fraction_available_mem 0.9

    for file in \$(ls ${sample_id})
    do
    [ -f ${sample_id}/"\$file" ] && mv ${sample_id}/"\$file" ${sample_id}_\${file}
    done

    echo "Finished PSV_TRIMQC" `date +%H-%M-%S`
    """
}

// Process - Align reads
process PSV_ALIGN {
    publishDir "${params.output_path}/alignments", mode: "copy"
    conda params.conda.psv_env

    executor "slurm"
    clusterOptions "-q gp_bscls -c ${params.psv_align.threads}"
    time params.psv_align.time
    array params.psv_align.array_size
    cpus params.psv_align.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(reference_file), \
        path(index), \
        path(reads_fwd), \
        path(reads_rev)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}_sorted.bam"), \
        path("${sample_id}_sorted.bam.bai")

    script:
    """
    echo "Starting PSV_ALIGN" `date +%H-%M-%S`
    ${params.psv_exe} align_reads \
        -o ${sample_id} \
        -f1 ${reads_fwd} \
        -f2 ${reads_rev} \
        -r ${reference_file} \
        -thr ${task.cpus} \
        --skip_marking_duplicates \
        --fraction_available_mem 0.9
    
    for file in \$(ls ${sample_id})
    do
    [ -f ${sample_id}/"\$file" ] && mv ${sample_id}/"\$file" ${sample_id}_\${file}
    done

    mv ${sample_id}_aligned_reads.bam.sorted ${sample_id}_sorted.bam
    mv ${sample_id}_aligned_reads.bam.sorted.bai ${sample_id}_sorted.bam.bai

    echo "Finished PSV_ALIGN" `date +%H-%M-%S`
    """
}

// Process - Subset alignment
process SUBSET_ALINGMENT{
    publishDir "${params.output_path}/alignment_subsets", mode: "copy"

    executor "local"
    time { 10.min * (task.attempt ** 2) }
    errorStrategy "retry"
    cpus 14
    maxRetries 2

    input:
    tuple \
        val(sample_id), \
        path(alignment), \
        path(index), \
        path(subset_bed)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}_sorted_subset_${subset_bed.simpleName}.bam"), \
        path("${sample_id}_sorted_subset_${subset_bed.simpleName}.bam.bai")

    script:
    """
    echo "Starting SUBSET_ALINGMENT" `date +%H-%M-%S`
    # Subset
    bedtools intersect -a ${alignment} -b ${subset_bed} > ${sample_id}_sorted_subset_${subset_bed.simpleName}.bam

    # Index
    samtools index ${sample_id}_sorted_subset_${subset_bed.simpleName}.bam ${sample_id}_sorted_subset_${subset_bed.simpleName}.bam.bai

    echo "Finished SUBSET_ALINGMENT" `date +%H-%M-%S`
    """
}

// Process - Alignment QC
process BAMSTATS {
    publishDir "${params.output_path}/bamstats", mode: "copy"
    conda params.conda.kat_env
    
    input:
    tuple \
        val(sample_id), \
        path(alignment), \
        path(alignment_index)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}_sam_stats.txt"), \
        path("${sample_id}_cov_stats.txt")

    script:
    """
    echo "Starting BAMSTATS" `date +%H-%M-%S`
    samtools stats -@ ${params.threads.general} ${alignment} > ${sample_id}_sam_stats.txt
    samtools coverage ${alignment} > ${sample_id}_cov_stats.txt
    echo "Finished BAMSTATS" `date +%H-%M-%S`
    """
}