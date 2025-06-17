// Process - Reference indexing with everything
process INDEX_FULL {
    executor "local"
    publishDir "${params.data_path}/references", mode: "copy"

    input:
    tuple \
        val(sample_id), \
        path(reference_file)

    output:
    tuple \
        val(sample_id), \
        path(reference_file), \
        path("${reference_file}.{amb,ann,bwt,pac,sa,fai,dict,bed}")

    script:
    """
    echo "Starting INDEX_FULL" `date +%H-%M-%S`
    bwa index $reference_file
    samtools faidx $reference_file > ${reference_file}.fai
    samtools dict $reference_file > ${reference_file}.dict
    awk -F'\\t' '{printf "%s\\t1\\t%s\\t%s\\n", \$1, \$2, \$1}' \
        ${reference_file}.fai > \
        ${reference_file}.bed
    echo "Finished INDEX_FULL" `date +%H-%M-%S`
    """
}

// Process - Reference indexing with everything
process PSV_REPEATS {
    publishDir "${params.data_path}/references", mode: "copy"
    conda params.conda.psv_env

    executor "slurm"
    clusterOptions "-q ${params.queue} -c ${params.psv_repeats.threads}"
    time params.psv_repeats.time
    array params.psv_repeats.array_size
    cpus params.psv_repeats.threads
    errorStrategy "ignore"
    input:
    tuple \
        val(sample_id), \
        path(reference_file)

    output:
    tuple \
        val(sample_id), \
        path("${reference_file.simpleName}_combined_repeats.tab")

    script:
    """
    echo "Starting PSV_REPEATS" `date +%H-%M-%S`
    ${params.psv_exe} infer_repeats \
        -o ${reference_file.simpleName} \
        -r ${reference_file} \
        -thr ${task.cpus} \
        --fraction_available_mem 0.9
    
    for file in \$(ls ${reference_file.simpleName})
    do
    [ -f ${reference_file.simpleName}/"\$file" ] && mv ${reference_file.simpleName}/"\$file" ${reference_file.simpleName}_\${file}
    done
    echo "Finished PSV_REPEATS" `date +%H-%M-%S`
    """
}