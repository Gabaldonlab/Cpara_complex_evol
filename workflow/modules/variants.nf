// Process - PerSVade call_small_vars
process PSV_CALLSNPS {
    publishDir "${params.output_path}/variant_calling", mode: "copy"
    conda params.conda.psv_env

    executor "slurm"
    clusterOptions "-q gp_bscls -c ${params.psv_callsnps.threads}"
    time params.psv_callsnps.time
    array params.psv_callsnps.array_size
    cpus params.psv_callsnps.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(reference_file), \
        path(index), \
        path(repeats_file), \
        path(reads_fwd), \
        path(reads_rev), \
        path(alignment), \
        path(alignment_index)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}_merged_vcfs_allVars_ploidy2.vcf"), \
        path("${sample_id}_variants_atLeast1PASS_ploidy2.vcf"), \
        path("${sample_id}_variants_atLeast2PASS_ploidy2.vcf"), \
        path("${sample_id}_variants_atLeast3PASS_ploidy2.vcf")

    script:
    """
    echo "Starting PSV_CALLSNPS" `date +%H-%M-%S`
    module unload bcftools
    python ${params.psv_exe} \
        call_small_variants \
        -o ${sample_id} \
        -r ${reference_file} \
        -sbam ${alignment} \
        --repeats_file ${repeats_file} \
        --callers bcftools,freebayes,HaplotypeCaller \
        --threads ${task.cpus} -p 2 \
        --min_AF 0.25 \
        --min_coverage 30 \
        --fraction_available_mem 0.95
   
    for file in \$(ls ${sample_id})
    do
    [ -f ${sample_id}/"\$file" ] && mv ${sample_id}/"\$file" ${sample_id}_\${file}
    done

    rm *withMultiAlt.vcf

    echo "Finished PSV_CALLSNPS" `date +%H-%M-%S`
    """
}

// Process - Compress SNPs
process COMPRESS_VCF {
    publishDir "${params.output_path}/variant_calling", mode: "copy"

    executor "local"
    time params.compress_vcf.time
    cpus params.compress_vcf.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(vcf)

    output:
    tuple \
        val(sample_id), \
        path("${vcf}.gz"), \
        path("${vcf}.gz.csi")

    script:
    """
    echo "Starting COMPRESS_VCF" `date +%H-%M-%S`
    # Compress
    bgzip -c -@ ${task.cpus} --keep ${vcf} > ${vcf}.gz
    # Index
    bcftools index -f ${vcf}.gz
    echo "Finished COMPRESS_VCF" `date +%H-%M-%S`
    """    
}

// Process - Compress SNPs
process DECOMPRESS_VCF {
    executor "local"
    time params.compress_vcf.time
    cpus params.compress_vcf.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(vcf_gz)

    output:
    tuple \
        val(sample_id), \
        path("*.vcf")

    script:
    """
    echo "Starting DECOMPRESS_VCF" `date +%H-%M-%S`
    # Decompress
    bgzip -d -@ ${task.cpus} ${vcf_gz}
    echo "Finished DECOMPRESS_VCF" `date +%H-%M-%S`
    """    
}

// Process - Subset background SNPs
process SUBSET_SNPS {
    publishDir "${params.output_path}/variant_subset", mode: "copy"

    executor "local"
    //clusterOptions "-q gp_bscls -c ${params.subset_snps.threads}"
    time params.subset_snps.time
    cpus params.subset_snps.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(vcf_gz), \
        path(vcf_csi), \
        val(background_id), \
        path(vcf_gz_background), \
        path(vcf_csi_background)

    output:
    tuple \
        val(sample_id), \
        path("${vcf_gz.simpleName}_wo_${background_id}.vcf.gz"), \
        path("${vcf_gz.simpleName}_wo_${background_id}.vcf.gz.csi")

    script:
    prefix = "${vcf_gz.simpleName}_wo_${background_id}"
    """
    echo "Starting SUBSET_SNPS" `date +%H-%M-%S`

    # Intersect
    bcftools isec -p ./ ${vcf_gz} ${vcf_gz_background}
    mv 0000.vcf ${prefix}.vcf
    # Compress
    bgzip -c -@ ${task.cpus} --keep \
        ${prefix}.vcf > \
        ${prefix}.vcf.gz
    # Index
    bcftools index -f ${prefix}.vcf.gz
    # Clean
    rm 000*.vcf ${prefix}.vcf

    echo "Finished SUBSET_SNPS" `date +%H-%M-%S`
    """     
}

// Process - Split variants by genotype
process SPLIT_VARIANTS{
    publishDir "${params.output_path}/variant_loh", mode: "copy"

    executor "local"
    //clusterOptions "-q gp_bscls -c ${params.subset_snps.threads}"
    time params.subset_snps.time
    cpus params.subset_snps.threads
    errorStrategy "ignore"

    input:
    tuple \
        val(sample_id), \
        path(vcf)

    output:
    tuple \
        val(sample_id), \
        path("${vcf.simpleName}_homVars.vcf.gz"), \
        path("${vcf.simpleName}_homVars.vcf.gz.csi"), emit: homVars
    tuple \
        val(sample_id), \
        path("${vcf.simpleName}_hetVars.vcf.gz"), \
        path("${vcf.simpleName}_hetVars.vcf.gz.csi"), emit: hetVars

    script:
    """
    echo "Starting SPLIT_VARIANTS" `date +%H-%M-%S`

    bcftools view -i 'GT="hom"' ${vcf} \
        -o ${vcf.simpleName}_homVars.vcf
    bcftools view -i 'GT="het"' ${vcf} \
        -o ${vcf.simpleName}_hetVars.vcf
    
    # Compress
    bgzip -c -@ ${task.cpus} \
        ${vcf.simpleName}_homVars.vcf > \
        ${vcf.simpleName}_homVars.vcf.gz
    bgzip -c -@ ${task.cpus} \
        ${vcf.simpleName}_hetVars.vcf > \
        ${vcf.simpleName}_hetVars.vcf.gz
    
    # Index
    bcftools index -f ${vcf.simpleName}_homVars.vcf.gz
    bcftools index -f ${vcf.simpleName}_hetVars.vcf.gz

    echo "Finished SPLIT_VARIANTS" `date +%H-%M-%S`
    """       
}

// Process - Detect LOH variants
process LOH_VARIANTS{
    publishDir "${params.output_path}/variant_loh", mode: "copy"
    
    executor "slurm"
    errorStrategy "ignore"
    clusterOptions "-q gp_bscls -c ${params.subset_snps.threads}"
    time params.subset_snps.time
    array params.jloh_extract.array_size
    cpus params.subset_snps.threads

    input:
    tuple \
        val(sample_id), \
        path(vcf_gz_all), \
        path(vcf_csi_all), \
        path(vcf_gz_homo), \
        path(vcf_csi_homo), \
        path(vcf_gz_heter), \
        path(vcf_csi_heter), \
        path(vcf_gz_bg_homo), \
        path(vcf_csi_bg_homo), \
        path(vcf_gz_bg_heter), \
        path(vcf_csi_bg_heter), \
        path(vcf_gz_t0_heter), \
        path(vcf_csi_t0heter)

    output:
    tuple \
        val(sample_id), \
        path("${sample_id}_lohAltVars.vcf"), \
        path("${sample_id}_lohRefVars.vcf"), \
        path("${sample_id}_nonLohVars.vcf")

    script:
    
    """
    echo "Starting LOH_VARIANTS" `date +%H-%M-%S`

    # Remove homozygous T0

    bcftools isec -p ./ ${vcf_gz_homo} ${vcf_gz_bg_homo}
    mv 0000.vcf ${sample_id}_homVars_wo_t0.vcf
    bcftools isec -p ./ ${vcf_gz_heter} ${vcf_gz_bg_homo}
    mv 0000.vcf ${sample_id}_hetVars_wo_t0.vcf

    # Compress
    bgzip -c -@ ${task.cpus} \
        ${sample_id}_homVars_wo_t0.vcf > \
        ${sample_id}_homVars_wo_t0.vcf.gz
    bgzip -c -@ ${task.cpus} \
        ${sample_id}_hetVars_wo_t0.vcf > \
        ${sample_id}_hetVars_wo_t0.vcf.gz

    # Index
    bcftools index -f ${sample_id}_homVars_wo_t0.vcf.gz
    bcftools index -f ${sample_id}_hetVars_wo_t0.vcf.gz

    # Variants that had LOH to ALT
    # 3PASS hetero at T0 and 3PASS homo in the evolved
    bcftools isec -p ./ ${sample_id}_homVars_wo_t0.vcf.gz ${vcf_gz_t0_heter}
    mv 0002.vcf ${sample_id}_lohAltVars.vcf

    # Variants that had LOH to REF
    # 3PASS hetero at T0 but no SNP in the allMerged evolved
    bcftools isec -p ./ ${vcf_gz_t0_heter} ${vcf_gz_all} 
    mv 0000.vcf ${sample_id}_lohRefVars.vcf

    # Variants that did not have LOH
    # 3PASS hetero at T0 and 3PASS hetero in the evolved
    bcftools isec -p ./ ${sample_id}_hetVars_wo_t0.vcf.gz ${vcf_gz_t0_heter}
    mv 0002.vcf ${sample_id}_nonLohVars.vcf

    # Remove decimals in AF
    sed -i \
        's/\\.6667:/:/g; s/\\.667:/:/g; s/\\.3333:/:/g; s/\\.333:/:/g; s/\\.33:/:/g; s/\\.3:/:/g' \
        ${sample_id}_lohAltVars.vcf
    sed -i \
        's/\\.6667:/:/g; s/\\.667:/:/g; s/\\.3333:/:/g; s/\\.333:/:/g; s/\\.33:/:/g; s/\\.3:/:/g' \
        ${sample_id}_lohRefVars.vcf
    sed -i \
        's/\\.6667:/:/g; s/\\.667:/:/g; s/\\.3333:/:/g; s/\\.333:/:/g; s/\\.33:/:/g; s/\\.3:/:/g' \
        ${sample_id}_nonLohVars.vcf

    echo "Finished LOH_VARIANTS" `date +%H-%M-%S`
    """       
}