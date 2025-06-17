// Process - Align reads
process DEPOOL_READS {
    publishDir "${params.data_path}/depooled_reads", mode: "copy"
    conda params.conda.depool_env

    executor "slurm"
    clusterOptions "-q gp_bscls -c ${params.threads.depool}"
    time { 30.min * task.attempt }
    errorStrategy "retry"
    array params.array_sizes.depool
    cpus params.threads.depool
    maxRetries 2

    input:
    tuple \
        val(pool), \
        val(sample_id), \
        val(pooled_genomes), \
        val(pooled_samples), \
        path(pooled_fastq1), \
        path(pooled_fastq2)

    output:
    tuple \
        val(sample_id), \
        path("${pool}/${sample_id}.1.fastq.gz"), \
        path("${pool}/${sample_id}.2.fastq.gz"), \
        path("${pool}/${pool}_read_splitting_statistics.tbl")

    script:
    """
    echo "Starting DEPOOL_READS" `date +%H-%M-%S`
    # Genomes with mitochondria (candu NO)
    CANOR=${params.depool_genomes.CANOR}
    CANME=${params.depool_genomes.CANME}
    CANGL=${params.depool_genomes.CANGL}
    CANAL=${params.depool_genomes.CANAL}
    CANAU=${params.depool_genomes.CANAU}
    CANDU=${params.depool_genomes.CANDU}
    SACCE=${params.depool_genomes.SACCE}
    SACUV=${params.depool_genomes.SACUV}
    SACKR=${params.depool_genomes.SACKR}
    SACCExSACKR=${params.depool_genomes.SACCExSACKR}
    SACCExSACUV=${params.depool_genomes.SACCExSACUV}
    SACCExSACAR=${params.depool_genomes.SACCExSACAR}

    python ${params.depool_script} \
        --fastq1 ${pooled_fastq1} \
        --fastq2 ${pooled_fastq2} \
        --outdir ${pool} \
        --threads ${params.threads.depool} \
        --refGenomes ${pooled_genomes} \
        --sampleNames ${pooled_samples}

    rm ${pool}/tmp_*
    mv ${pool}/read_splitting_statistics.tbl ${pool}/${pool}_read_splitting_statistics.tbl

    echo "Finished DEPOOL_READS" `date +%H-%Mq-%S`
    """
}
