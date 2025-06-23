// Process - PerSVade call_CNVs
process PSV_CALLCNVS {
    conda params.conda.psv_env_miki

    executor "slurm"
    clusterOptions "-q ${params.queue} -c ${params.psv_callcnvs.threads}"
    time params.psv_callcnvs.time
    array params.psv_callcnvs.array_size
    cpus params.psv_callcnvs.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(reference_file), \
        path(index), \
        path(repeats_file), \
        val(mitochondria), \
        path(reads_fwd), \
        path(reads_rev), \
        path(alignment), \
        path(alignment_index)

    output:
    tuple \
        val(sample_id), \
        path(sample_id)

    script:
    """
    echo "Starting PSV_CALLCNVS" `date +%H-%M-%S`
    python ${params.psv_exe_miki} \
        call_CNVs \
        -o ${sample_id} \
        -r ${reference_file} \
        -sbam ${alignment} \
        --mitochondrial_chromosome ${mitochondria} \
        --repeats_file ${repeats_file} \
        --threads ${task.cpus} -p 2 \
        --cnv_calling_algs ${params.psv_callcnvs.algs} \
        --window_size_CNVcalling ${params.psv_callcnvs.window_size} \
        --fraction_available_mem 0.95

    echo "Finished PSV_CALLCNVS" `date +%H-%M-%S`
    """
}

// Process - Filter CNV calling based on coverage
process CNV_FILTER {
    executor "local"
    
    input:
    tuple \
        val(sample_id), \
        path(cnv_outdir)

    output:
    tuple \
        val(sample_id), \
        path(sample_id)

    script:
    """
    echo "Starting CNV_FILTER" `date +%H-%M-%S`

    mv ${sample_id}/final_CNVcalling.tab ${sample_id}/unfiltered_final_CNVcalling.tab
    Rscript ${projectDir}/templates/cnv_covfilt.r \
        ${params.cnv_filter.max_del_cov} ${params.cnv_filter.min_dup_cov}
    mv final_CNVcalling.tab ${sample_id}/

    echo "Finished CNV_FILTER" `date +%H-%M-%S`
    """
}

// Process - PerSVade optimize_parameters
process PSV_OPTIMIZE {
    publishDir "${params.output_path}/sv_optimize", mode: "copy"
    conda params.conda.psv_env_miki

    executor "slurm"
    clusterOptions "-q ${params.queue}"
    time params.psv_optimize.time
    array params.psv_optimize.array_size
    cpus params.psv_optimize.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(reference_file), \
        path(index), \
        path(repeats_file), \
        val(mitochondria), \
        path(alignment), \
        path(alignment_index)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}_optimized_parameters.json"), \
        path("${sample_id}")

    script:
    """
    echo "Starting PSV_OPTIMIZE" `date +%H-%M-%S`
    python ${params.psv_exe_miki} \
        optimize_parameters \
        -o ${sample_id} \
        -r ${reference_file} \
        -sbam ${alignment} \
        --repeats_file ${repeats_file} \
        --mitochondrial_chromosome ${mitochondria} \
        --threads ${task.cpus} \
        --regions_SVsimulations ${params.psv_optimize.regions} \
        --simulation_ploidies ${params.psv_optimize.ploidies} \
        --keep_simulation_files \
        --fraction_available_mem 0.95
   
    for file in \$(ls ${sample_id})
    do
    [ -f ${sample_id}/"\$file" ] && mv ${sample_id}/"\$file" ${sample_id}_\${file}
    done

    echo "Finished PSV_OPTIMIZE" `date +%H-%M-%S`
    """
}

// Process - PerSVade analyze_SV_parameters
process PSV_SVPARAMS {
    publishDir "${params.output_path}/sv_params", mode: "copy"
    conda params.conda.psv_env_miki

    executor "slurm"
    clusterOptions "-q ${params.queue} -c ${params.psv_svparams.threads}"
    time params.psv_svparams.time
    array params.psv_svparams.array_size
    cpus params.psv_svparams.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(reference_id), \
        path(reference_file), \
        path(index), \
        path(repeats_file), \
        val(mitochondria), \
        path(opt_jsons), \
        path(opt_dirs)

    output:
    path(reference_id)

    script:
    """
    echo "Starting PSV_SVPARAMS" `date +%H-%M-%S`
    module unload bcftools

    python ${params.psv_exe_miki} \
        analyze_SV_parameters \
        -o ${reference_id} \
        --SV_parameters ${opt_jsons} \
        --optimize_parameters_outdirs ${opt_dirs} \
        -r ${reference_file} \
        -mchr ${mitochondria} \
        --repeats_file ${repeats_file} \
        -thr ${task.cpus} \
        --fraction_available_mem 0.95

    echo "Finished PSV_SVPARAMS" `date +%H-%M-%S`
    """
}

// Process - PerSVade call_SVs
process PSV_CALLSVS {
    publishDir "${params.output_path}/sv_calling", mode: "copy"
    conda params.conda.psv_env

    executor "slurm"
    clusterOptions "-q ${params.queue} -c ${params.psv_callsvs.threads}"
    time params.psv_callsvs.time
    array params.psv_callsvs.array_size
    cpus params.psv_callsvs.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(reference_file), \
        path(index), \
        path(repeats_file), \
        val(mitochondria), \
        path(opt_params), \
        path(alignment), \
        path(alignment_index)

    output:
    tuple \
        val(sample_id), \
        path(sample_id)

    script:
    """
    echo "Starting PSV_CALLSVS" `date +%H-%M-%S`
    module unload bcftools

    python ${params.psv_exe} \
        call_SVs \
        -o ${sample_id} \
        -r ${reference_file} \
        -sbam ${alignment} \
        --repeats_file ${repeats_file} \
        --mitochondrial_chromosome ${mitochondria} \
        --threads ${task.cpus} \
        --SVcalling_parameters ${opt_params} \
        --fraction_available_mem 0.95

    echo "Finished PSV_CALLSVS" `date +%H-%M-%S`
    """
}

// Process - PerSVade integrate_CNVs_SVs
process PSV_INTEGRATE {
   publishDir "${params.output_path}/sv_integrated", mode: "copy"
    conda params.conda.psv_env

    executor "slurm"
    clusterOptions "-q ${params.queue} -c ${params.psv_integrate.threads}"
    time { params.psv_integrate.time * task.attempt }
    array params.psv_integrate.array_size
    cpus params.psv_integrate.threads
    errorStrategy "retry"
    maxRetries 2

    input:
    tuple \
        val(sample_id), \
        path(reference_file), \
        path(index), \
        path(repeats_file), \
        val(mitochondria), \
        path(gff), \
        path(bam), \
        path(bai), \
        path(sv_output, stageAs: "svs"), \
        path(cnv_output, stageAs: "cnvs")

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}")

    script:
    """
    echo "Starting PSV_INTEGRATE" `date +%H-%M-%S`
    python ${params.psv_exe} \
        integrate_SV_CNV_calls \
        -o ${sample_id} \
        -r ${reference_file} \
        --repeats_file ${repeats_file} \
        --mitochondrial_chromosome ${mitochondria} \
        -sbam ${bam} -p 2 \
        --outdir_callSVs ${sv_output} \
        --outdir_callCNVs ${cnv_output} \
        --threads ${task.cpus} \
        --fraction_available_mem 0.95
    
    echo "Finished PSV_INTEGRATE" `date +%H-%M-%S`
    """
}

// Process - PerSVade integrate_several_samples
process PSV_COMPARE {
    publishDir "${params.output_path}/sv_background", mode: "copy"
    conda params.conda.psv_env_miki

    executor "slurm"
    clusterOptions "-q ${params.queue} -c ${params.psv_compare.threads}"
    time { params.psv_compare.time * task.attempt }
    array params.psv_compare.array_size
    cpus params.psv_compare.threads
    errorStrategy "retry"
    maxRetries 3

    input:
    tuple \
        val(reference_id), \
        path(reference_file), \
        path(index), \
        path(repeats_file), \
        val(mitochondria), \
        path(sv_tables), \
        path(comp_table)

    output:
    path(reference_id)

    script:
    """
    echo "Starting PSV_COMPARE" `date +%H-%M-%S`
    
    python ${params.psv_exe_miki} \
        integrate_several_samples \
        -o ${reference_id} \
        --paths_table ${sv_tables} \
        --table_comparisons ${comp_table} \
        -r ${reference_file} \
        -mchr ${mitochondria} \
        --repeats_file ${repeats_file} \
        -p 2 \
        --tol_bp ${params.psv_compare.tolbp} \
        -thr ${task.cpus} \
        --fraction_available_mem 0.95

    echo "Finished PSV_COMPARE" `date +%H-%M-%S`
    """
}