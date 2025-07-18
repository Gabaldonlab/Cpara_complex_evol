// Log info
log.info """\
        ==========================
        LOG INFO
        ==========================
        output path  : ${params.output_path}
        --
        run as       : ${workflow.commandLine}
        started at   : ${workflow.start}
        config files : ${workflow.configFiles}
        container    : ${workflow.containerEngine}:${workflow.container}
        ==========================
         """
        .stripIndent()

// Modules
include { DEPOOL_READS } from './modules/depool.nf'
include { INDEX_FULL; PSV_REPEATS } from './modules/references.nf'
include { 
    PSV_CALLSNPS; COMPRESS_VCF; DECOMPRESS_VCF; DECOMPRESS_VCF as DECOMPRESS_NEW;
    DECOMPRESS_VCF as DECOMPRESS_NEWHOM; SUBSET_SNPS as SUBSET_SNPS_T0;
    SUBSET_SNPS as SUBSET_SNPS_YPD; SPLIT_VARIANTS as SPLIT_ALLVARS;
    SPLIT_VARIANTS as SPLIT_3PASS; SPLIT_VARIANTS as SPLIT_WOBG;
    SPLIT_VARIANTS as SPLIT_NEW;
    LOH_VARIANTS
} from './modules/variants.nf'
include {
    PSV_TRIMQC; PSV_ALIGN; SUBSET_ALINGMENT; BAMSTATS
} from './modules/reads.nf'
include {
    JLOH_STATS; PARSE_JLOHSTATS; JLOH_EXTRACT; JLOH_PLOT; MERGE_BEDS;
    BED_INTERSECT as BED_INTERSECT_T0; BED_INTERSECT as BED_INTERSECT_YPD;
    BED_VCF_INTERSECT; JLOH_PLOT_VARS; BED_LOH_STATS
} from './modules/jloh.nf'
include {
    PSV_CALLCNVS; PSV_OPTIMIZE; PSV_CALLSVS; PSV_SVPARAMS; PSV_INTEGRATE;
    PSV_COMPARE; CNV_FILTER; PSV_SVPARAMS as PSV_MANPARAMS
} from './modules/structural.nf'
include {
    PSV_ANNOTSNPS; PSV_ANNOTSVS; MERGE_GFFSNP; SNP_PLOT; MERGE_GFFLOH
} from './modules/annotate.nf'

// Channels of raw data
PooledReads = Channel.fromFilePairs(
    "${params.pools_path}/*_{1,2}.fastq.gz",
    flat: true
)
Metadata = Channel.fromPath("${params.metadata_path}/estrain_metadata_workflow.tsv")
    .splitCsv( skip: 1, sep: "\t", strip: true)
    .filter{ it[8] == "false" }

References = Channel.fromPath("${params.data_path}/references/*.{fna,fasta,fa}")
    .map{ [it.simpleName, it] }

ReferenceGff = Channel.fromPath("${params.data_path}/references/*.gff")
    .map{ [it.simpleName, it] }

Mitochondria = Channel.of(
    ["CANOR_CO90125_mit_MCO471", "DQ026513.2"],
    ["CANME_BP57_mit_BP57", "JQ062879.1"]
)
SubsetBed = Channel.of(
    ["CANOR_CO90125_mit_MCO471", "${params.data_path}/subset_CANOR_FKS1.bed"],
    ["CANOR_CO90125_mit_MCO471", "${params.data_path}/subset_CANOR_MRR1.bed"],
    ["CANME_BP57_mit_BP57", "${params.data_path}/subset_CANME_FKS1.bed"],
    ["CANME_BP57_mit_BP57", "${params.data_path}/subset_CANME_MRR1.bed"]
)

// Workflow to index and generate repeats table of reference genomes
workflow REFERENCES {
    main:
        INDEX_FULL(References)
        PSV_REPEATS(References)
}

// Channels of the output of REFERENCES
ReferenceRepeats = Channel.fromPath(
    "${params.data_path}/references/*combined_repeats.tab")
    .map{ [it.getName().replaceAll("_combined_repeats.tab", ""), it] }

ReferenceIndices = Channel.fromPath(
    "${params.data_path}/references/*.{amb,ann,bwt,pac,sa,fai,dict,bed}")
    .map{ [it.getName().replaceAll("\\..*", ""), it] }
    .groupTuple()

IndexedRefereces = References.join(ReferenceIndices)

// Workflow to depool reads
workflow DEPOOL {
    main:
        MatchSamplePool = Metadata
            .map{ [it[7], it[0], it[9], it[10]] }
            .join(PooledReads)
        DEPOOL_READS(MatchSamplePool)
}

// Channel of depooled reads
RawReads = Channel.fromFilePairs(
    "${params.reads_path}/*/EPOM*.{1,2}.fastq.gz",
    flat: true
)

// Workflow to align and do quality checks
workflow ALIGN {
    main:
        FilteredReads = PSV_TRIMQC(RawReads)
        MatchSampleRef = IndexedRefereces
            .join(ReferenceRepeats)
            .cross(Metadata.map{ [it[5], it[0]] })
            .map{ [it[1][1], it[0][1], it[0][2]] }
            .join(FilteredReads.map{ [it[0], it[1], it[2]] })
        Alignment = PSV_ALIGN(MatchSampleRef)
        BamStats = BAMSTATS(Alignment)
}

// Channels of the output of ALIGN
FilteredReads = Channel.fromFilePairs(
    "${params.output_path}/reads_filtered_and_qc/EPOM*{1,2}.fastq.gz",
    flat: true)
    .map{ [it[0].replaceAll("_trimmed_reads", ""), it[1], it[2]] }

Alignment_bams = Channel.fromPath(
    "${params.output_path}/alignments/*.bam"
    )
    .map{ [it.name.replaceAll("_sorted.bam", ""), it] }

Alignment_indices = Channel.fromPath(
    "${params.output_path}/alignments/*.bam.bai"
    )
    .map{ [it.name.replaceAll("_sorted.bam.bai", ""), it] }

Alignments = Alignment_bams.join(Alignment_indices)

// Workflow for small variant calling and background removal
workflow VARIANTS {
    main:
        MatchToCall = IndexedRefereces
            .join(ReferenceRepeats)
            .cross(Metadata.map{ [it[5], it[0]] })
            .map{ [it[1][1], it[0][1], it[0][2], it[0][3]] }
            .join(FilteredReads)
            .join(Alignments)
        Variants = PSV_CALLSNPS(MatchToCall)
        // Get the variants and the condition
        VariantsAndCondition = Variants.join(Metadata.map{ [it[0], it[1]] })
        // Select the all_merged for background
        T0Variants = VariantsAndCondition.filter{ it[5] == "T0" }
            .map{ [it[0], it[1]] }
        YpdVariants = VariantsAndCondition.filter{ it[5] == "Y" }
            .map{ [it[0], it[1]] }
       // Select the 3PASS for evolved
        EvolvedVariants = VariantsAndCondition
            .filter{ it[5] != "T0" && it[5] != "Y"}
            .map{ [it[0], it[4]] }
        // Compress files
        CompressedVcfs = COMPRESS_VCF(
            EvolvedVariants.mix(YpdVariants).mix(T0Variants)
        )
        // Match with clade
        VariantsAndClade = CompressedVcfs
            .combine(Metadata.map{ [it[0], it[1], it[4]] }, by: 0)
        // Separate T0, YPD and evolved
        T0CompressedVcfs = VariantsAndClade.filter{ it[3] == "T0" }
        YpdCompressedVcfs = VariantsAndClade.filter{ it[3] == "Y" }
        EvolvedCompressedVcfs = VariantsAndClade
            .filter{ it[3] != "T0" && it[3] != "Y"}
        // Remove background matching T0 by clade
        EvolvedMatchT0 = EvolvedCompressedVcfs
            .combine(T0CompressedVcfs, by: 4)
            .map{ [it[1], it[2], it[3], it[5], it[6], it[7]] }
        EvolvedWoT0 = SUBSET_SNPS_T0(EvolvedMatchT0)
        // Remove background matching YPD by clade
        EvolvedWoT0AndClade = EvolvedWoT0
            .combine(Metadata.map{ [it[0], it[1], it[4]] }, by: 0)
        EvolvedMatchYPD = EvolvedWoT0AndClade
            .combine(YpdCompressedVcfs, by: 4)
            .map{ [it[1], it[2], it[3], it[5], it[6], it[7]] }
        EvolvedWoYPD = SUBSET_SNPS_YPD(EvolvedMatchYPD)
}

// Channels from variants
VariantsAllVars = Channel.fromPath(
    "${params.output_path}/variant_calling/*allVars*vcf"
    )
    .map{ [it.name.split("_")[0..2].join("_"), it] }

Variants2Pass = Channel.fromPath(
    "${params.output_path}/variant_calling/*2PASS*vcf"
    )
    .map{ [it.name.split("_")[0..2].join("_"), it] }

Variants3Pass = Channel.fromPath(
    "${params.output_path}/variant_calling/*3PASS*vcf"
    )
    .map{ [it.name.split("_")[0..2].join("_"), it] }

VariantsNew = Channel.fromPath(
    "${params.output_path}/variant_subset/*wo*Y*gz"
    )
    .map{ [it.name.split("_")[0..2].join("_"), it] }

// Workflow for JLOH and LOH analyses
workflow LOH {
    main:
        JlohStats = JLOH_STATS(Variants2Pass)
        JlohStatsParsed = PARSE_JLOHSTATS(JlohStats)
        JlohExtract = JLOH_EXTRACT( IndexedRefereces
            .cross(Metadata.map{ [it[5], it[0]] })
            .map{ [it[1][1], it[0][1], it[0][2]] }
            .join(Alignments)
            .join(Variants2Pass)
            .join(JlohStatsParsed) )
        // Merge uncalled and LOH regions used as background
        LohWithUnc = MERGE_BEDS(JlohExtract.map{ [it[0], it[1], it[4]] })
        // Match with clade and split the T0 and YPD
        UncWithClade = LohWithUnc
            .combine(Metadata.map{ [it[0], it[1], it[4]] }, by: 0)
        T0JlohUnc = UncWithClade.filter{ it[2] == "T0" }
        YpdJlohUnc = UncWithClade.filter{ it[2] == "Y" }
        // Separate T0 for plotting later
        LohT0 = JlohExtract
            .map{ [it[0], it[1]] }
            .combine(Metadata.map{ [it[0], it[1]] }, by: 0)
            .filter{ it[2] == "T0"}
            .map{ [it[0], it[1]] }
        // Separate the evolved and YPD strains
        EvolvedJlohBeds = JlohExtract
            .map{ [it[0], it[1]] }
            .join(LohWithUnc)
            .combine(Metadata.map{ [it[0], it[1], it[4]] }, by: 0)
            .filter{ it[3] != "T0"}
            .map{ [it[0], it[1], it[2], it[4], it[3]] }
        // Remove background matching T0 by clade
        EvolvedWoT0 = BED_INTERSECT_T0(EvolvedJlohBeds
            .combine(T0JlohUnc, by: 3)
            .map{ [it[1], it[2], it[3], it[5], it[6]] })
        // Get statistics for the LOH blocks
        BED_LOH_STATS(EvolvedWoT0
            .combine(Metadata.map{ [it[0], it[5]] }, by: 0)
            .map{ [it[2], it[0], it[1]] }
            .combine(ReferenceIndices.map{ [it[0], it[1][6]] }, by: 0)
            .map{ [it[1], it[2], it[3]] }
            )
        // Remove background matching YPD by clade
        // EvolvedWoYPD = BED_INTERSECT_YPD(EvolvedWoT0
        //     .join(LohWithUnc)
        //     .combine(Metadata.map{ [it[0], it[4]] }, by: 0)
        //     .combine(YpdJlohUnc, by: 3)
        //     .map{ [it[1], it[2], it[3], it[5], it[6]] })
        // Plot LOH
        JLOH_PLOT(JlohExtract.map{ [it[0], it[1], it[2], it[3]] }
            // Substitute all LOH with the LOH without background
            .join(EvolvedWoT0.concat(LohT0))
            .map{ [it[0], it[4], it[2], it[3]] }
            // Add reference and clade ids
            .combine(Metadata.map{ [it[0], it[5], it[4]] }, by: 0)
            .map{ [it[4], it[5], [it[1], it[2], it[3]]] }
            // Add reference index
            .combine(ReferenceIndices, by: 0)
            // Group by clade
            .groupTuple( by: 1 )
            .map { [it[1], it[2].flatten(), it[3][0]] }
        )
        // Detect variants that underwent LOH
        // Split 3PASS variants based on genotype
        SPLIT_3PASS(Variants3Pass)
        // Filter 3PASS variants
        Vars3PassHomEvol = SPLIT_3PASS.out.homVars
            .combine(Metadata.map{ [it[0], it[1], it[4]] }, by: 0)
            .filter{ it[3] != "T0" }
            .map{ [it[0], it[1], it[2]] }   
        Vars3PassHetEvol = SPLIT_3PASS.out.hetVars
            .combine(Metadata.map{ [it[0], it[1], it[4]] }, by: 0)
            .filter{ it[3] != "T0" }
            .map{ [it[0], it[1], it[2]] }   
        Vars3PassHetT0 = SPLIT_3PASS.out.hetVars
            .combine(Metadata.map{ [it[0], it[1], it[4]] }, by: 0)
            .filter{ it[3] == "T0" }
            .map{ [it[0], it[1], it[2]] }   
        // Split allMerged variants based on genotype only for T0
        SPLIT_ALLVARS(VariantsAllVars
            .combine(Metadata.map{ [it[0], it[1], it[4]] }, by: 0)
            .filter{ it[2] == "T0" }
            .map{ [it[0], it[1]] })
        // Compress the allVars variants for evolved strains
        VarsAllAllEvol = COMPRESS_VCF(VariantsAllVars
            .combine(Metadata.map{ [it[0], it[1], it[4]] }, by: 0)
            .filter{ it[2] != "T0" }
            .map{ [it[0], it[1]] })
        // Combine all variants from evolved strains and add clade
        MergeVarsEvol = VarsAllAllEvol
            .join(Vars3PassHomEvol)
            .join(Vars3PassHetEvol)
            .combine(Metadata.map{ [it[0], it[4]] }, by: 0)
            .map{ [it[7], it[0], it[1], it[2], it[3], it[4], it[5], it[6]] }
        // Combine all variants from T0 strains and add clade
        MergeVarsT0 = SPLIT_ALLVARS.out.homVars
            .join(SPLIT_ALLVARS.out.hetVars)
            .join(Vars3PassHetT0)
            .combine(Metadata.map{ [it[0], it[4]] }, by: 0)
            .map{ [it[7], it[0], it[1], it[2], it[3], it[4], it[5], it[6]] }
        // Associate each sample to its T0 by clade
        VariantsMatchLoh = MergeVarsEvol
            .combine(MergeVarsT0, by: 0)
            .map{ [it[1], it[2], it[3], it[4], it[5], it[6], it[7],
                it[9], it[10], it[11], it[12], it[13], it[14]] }
        // Get variants that underwent LOH towards ALT/REF or did not have LOH
        // taking into account the switches between genotype from T0 to evolved
        VariantsLoh = LOH_VARIANTS(VariantsMatchLoh)
        // Add the newly-acquired hom variants
        VariantsNewVcf = DECOMPRESS_NEW(VariantsNew)
        SPLIT_NEW(VariantsNewVcf)
        VarintsNewHom = DECOMPRESS_NEWHOM(SPLIT_NEW.out.homVars.map{ [it[0], it[1]] })
        JLOH_PLOT_VARS(JlohExtract.map{ [it[0], it[1], it[2], it[3]] }
            // Add reference and clade
            .combine(Metadata.map{ [it[0], it[5], it[4]] }, by: 0)
            // Add new variants and LOH variants for each sample
            .join(VariantsLoh)
            .join(VarintsNewHom)
            // Map as reference, clade, and beds and vcfs
            .map{ [it[4], it[5], [
                it[1], it[2], it[3], it[6], it[7], it[8], it[9]
            ]] }
            // Add reference index
            .combine(ReferenceIndices, by: 0)
            // Group by clade
            .groupTuple( by: 1 )
            // Map as clade, beds and vcfs, and one genome index
            .map { [it[1], it[2].flatten(), it[3][0]] }
        )
}

CnvCalling = Channel.fromPath(
    "${params.output_path}/cnv_calling/*", type: 'dir'
    )
    .map{ [it.name, it] }

// Workflow for structural variants
workflow STRUCTURAL {
    main:
        // This channel will be useful
        RefWithExtras = IndexedRefereces
            .join(ReferenceRepeats)
            .join(Mitochondria)
        // Do the CNV calling
        // CnvCalling = PSV_CALLCNVS(RefWithExtras
        //     .combine(Metadata.map{ [it[5], it[0]] }, by: 0)
        //     .map{ [it[5], it[1], it[2], it[3], it[4]] }
        //     .join(FilteredReads)
        //     .join(Alignments)
        // )
        CnvFiltered = CNV_FILTER(CnvCalling)
        // Optimise parameters for the SVs
        OptimizedParameters = PSV_OPTIMIZE(RefWithExtras
            .combine(Metadata.map{ [it[5], it[0]] }, by: 0)
            .map{ [it[5], it[1], it[2], it[3], it[4]] }
            .join(Alignments))
        // Save the files for the parameter analysis
        OptJsonCanor = OptimizedParameters
            .combine(Metadata.map{ [it[0], it[5]] }, by: 0)
            .filter{ it[3] == "CANOR_CO90125_mit_MCO471" }
            .map{ it -> "parameterID\tparameters_json\n${it[0]}\t${it[1]}\n" }
            .collectFile(
                name: 'optimized_parameters_CANOR.tsv',
                storeDir:"${params.output_path}/sv_params",
                keepHeader: true,
                skip: 1
            )
            .map{ ["CANOR_CO90125_mit_MCO471", it] }
        OptJsonCanme = OptimizedParameters
            .combine(Metadata.map{ [it[0], it[5]] }, by: 0)
            .filter{ it[3] == "CANME_BP57_mit_BP57" }
            .map{ it -> "parameterID\tparameters_json\n${it[0]}\t${it[1]}\n" }
            .collectFile(
                name: 'optimized_parameters_CANME.tsv',
                storeDir:"${params.output_path}/sv_params",
                keepHeader: true,
                skip: 1
            )
            .map{ ["CANME_BP57_mit_BP57", it] }
        OptDirCanor = OptimizedParameters
            .combine(Metadata.map{ [it[0], it[5]] }, by: 0)
            .filter{ it[3] == "CANOR_CO90125_mit_MCO471" }
            .map{ it -> "sampleID\toptimize_parameters_outdir\n${it[0]}\t${it[2]}\n" }
            .collectFile(
                name: 'optimization_output_CANOR.tsv',
                storeDir:"${params.output_path}/sv_params",
                keepHeader: true,
                skip: 1
            )        
            .map{ ["CANOR_CO90125_mit_MCO471", it] }
        OptDirCanme = OptimizedParameters
            .combine(Metadata.map{ [it[0], it[5]] }, by: 0)
            .filter{ it[3] == "CANME_BP57_mit_BP57" }
            .map{ it -> "sampleID\toptimize_parameters_outdir\n${it[0]}\t${it[2]}\n" }
            .collectFile(
                name: 'optimization_output_CANME.tsv',
                storeDir:"${params.output_path}/sv_params",
                keepHeader: true,
                skip: 1
            )
            .map{ ["CANME_BP57_mit_BP57", it] }
        
        
        // Create the channel for those files
        OptJson = OptJsonCanor.concat(OptJsonCanme)
        OptDir = OptDirCanor.concat(OptDirCanme)
        // Combine everything for the SV parameter analysis
        RefSvParams = RefWithExtras
            .join(OptJson)
            .join(OptDir)
        // Run the analysis of SV parameters,
        // after this step you'll have to manually choose the parameters
        // PSV_SVPARAMS(RefSvParams)
        // Create the channel for the manually chosen parameters
        OptParamsManual = Channel.of(
            ["CANOR_CO90125_mit_MCO471", "${params.output_path}/sv_params/opt_params_manual_CANOR.json"],
            ["CANME_BP57_mit_BP57", "${params.output_path}/sv_params/opt_params_manual_CANME.json"]
        )
        // Check the manually chosen parameters
        // These channels generate the files for the parameters, but they need to be later edited manually
        // ManJsonCanor = OptimizedParameters
        //     .combine(Metadata.map{ [it[0], it[5]] }, by: 0)
        //     .filter{ it[3] == "CANOR_CO90125_mit_MCO471" }
        //     .map{ it -> "parameterID\tparameters_json\n${it[0]}\t${params.output_path}/sv_params/opt_params_manual_CANOR.json\n" }
        //     .collectFile(
        //         name: 'manual_parameters_CANOR.tsv',
        //         storeDir:"${params.output_path}/sv_params",
        //         keepHeader: true,
        //         skip: 1
        //     )
        //     .map{ ["CANOR_CO90125_mit_MCO471", it] }
        // ManJsonCanme = OptimizedParameters
        //     .combine(Metadata.map{ [it[0], it[5]] }, by: 0)
        //     .filter{ it[3] == "CANME_BP57_mit_BP57" }
        //     .map{ it -> "parameterID\tparameters_json\n${it[0]}\t${params.output_path}/sv_params/opt_params_manual_CANME.json\n" }
        //     .collectFile(
        //         name: 'manual_parameters_CANME.tsv',
        //         storeDir:"${params.output_path}/sv_params",
        //         keepHeader: true,
        //         skip: 1
        //     )
        //     .map{ ["CANME_BP57_mit_BP57", it] }
        // // Create the channel for those files
        // ManJson = ManJsonCanor.concat(ManJsonCanme)
        // This channel contains the files above, but after a manual modification
        ManJson = Channel.of(
            ["CANOR_CO90125_mit_MCO471", "${params.output_path}/sv_params/manual_parameters_CANOR.csv"],
            ["CANME_BP57_mit_BP57", "${params.output_path}/sv_params/manual_parameters_CANME.csv"]
        )  
        // Combine everything for the SV parameter analysis
        ManSvParams = RefWithExtras
            .join(ManJson)
            .join(OptDir)
        // Run the analysis of the manual SV parameters
        PSV_MANPARAMS(ManSvParams)
        // Combine the references with the manually chosen parameters
        // and join with the samples to do the SV calling
        StrucuturalVariants = PSV_CALLSVS(RefWithExtras
            .join(OptParamsManual)
            .combine(Metadata.map{ [it[5], it[0]] }, by: 0)
            .map{ [it[6], it[1], it[2], it[3], it[4], it[5]] }
            .join(Alignments)
        )
        // Integrate CNVs and SVs
        IntegratedSvs = PSV_INTEGRATE(RefWithExtras
            .join(ReferenceGff)
            .combine(Metadata.map{ [it[5], it[0]] }, by: 0)
            .map{ [it[6], it[1], it[2], it[3], it[4], it[5]] }
            .join(Alignments)
            .join(StrucuturalVariants)
            .join(CnvFiltered)
        )
        // Create the channel for the SV path tables
        SvPathCanme = IntegratedSvs
            .join(StrucuturalVariants)
            .join(Alignments)
            .combine(Metadata.map{ [it[0], it[5]] }, by: 0)
            .filter{ it[5] == "CANME_BP57_mit_BP57" }
            .map{ it -> "sampleID\tintegrate_SV_CNV_calls_outdir\tcall_SVs_outdir\tsorted_bam\n${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}\n" }
            .collectFile(
                name: 'sv_paths_table_CANME.tsv',
                storeDir:"${params.output_path}/sv_background",
                keepHeader: true,
                skip: 1
            )
            .map{ ["CANME_BP57_mit_BP57", it] }
        SvPathCanor = IntegratedSvs
            .join(StrucuturalVariants)
            .join(Alignments)
            .combine(Metadata.map{ [it[0], it[5]] }, by: 0)
            .filter{ it[5] == "CANOR_CO90125_mit_MCO471" }
            .map{ it -> "sampleID\tintegrate_SV_CNV_calls_outdir\tcall_SVs_outdir\tsorted_bam\n${it[0]}\t${it[1]}\t${it[2]}\t${it[3]}\n" }
            .collectFile(
                name: 'sv_paths_table_CANOR.tsv',
                storeDir:"${params.output_path}/sv_background",
                keepHeader: true,
                skip: 1
            )
            .map{ ["CANOR_CO90125_mit_MCO471", it] }
        // Join channels
        SvPathTables = SvPathCanme.concat(SvPathCanor)
        // Create the channel for the comparison tables
        SvCompTables = Channel.of(
            ["CANOR_CO90125_mit_MCO471", "${params.output_path}/sv_background/table_comparisons_canor.csv"],
            ["CANME_BP57_mit_BP57", "${params.output_path}/sv_background/table_comparisons_canme.csv"]
        )       
        // Remove background from SVs
        SvWoBg = PSV_COMPARE(RefWithExtras
            .join(SvPathTables)
            .join(SvCompTables)
        )
}

SvOutDirs = Channel.fromPath(
    "${params.output_path}/sv_calling/*", type: 'dir'
    )
    .map{ [it.name, it] }

CnvOutDirs = Channel.fromPath(
    "${params.output_path}/cnv_calling/*", type: 'dir'
    )
    .map{ [it.name, it] }
    
IntegratedSvs = Channel.fromPath(
    "${params.output_path}/sv_integrated/*/SV_and_CNV_variant_calling.vcf",
    )
    .map{ [it.name.split("/")[0], it] }

LohBlocks = Channel.fromPath(
    "${params.output_path}/jloh_extract/EPOM*wo_EPOM_T0_{5,6,7,8,9,10,11,12}.bed"
    )
    .map{ [it.name.split("_")[0..2].join("_"), it] }


// Workflow to annotate all the called variants
workflow ANNOTATE {
    main:
        // Generate GFF?
        // This channel will be useful
        RefWithExtras = IndexedRefereces
            .join(ReferenceGff)
            .join(Mitochondria)
        VariantsNewVcf = DECOMPRESS_VCF(VariantsNew)
        // Annotate new SNPs with T0 and YPD removed
        AnnotSnps = PSV_ANNOTSNPS(RefWithExtras
            .combine(Metadata.map{ [it[5], it[0]] }, by: 0)
            .map{ [it[5], it[1], it[2], it[3], it[4]] }
            .join(VariantsNewVcf)
        )
        // Merge PSV annotations with the GFF
        // GffSnps = MERGE_GFFSNP(AnnotSnps
        //     .combine(Metadata.map{ [it[0], it[5]] }, by: 0)
        //     .map{ [it[2], it[1]] }
        //     .groupTuple(by: 0)
        //     .map{ [it[0], it[1].flatten()] }
        //     .join(ReferenceGff)
        // )
        // Merge acquired SNPs with LOH recombination blocks
        BED_VCF_INTERSECT(VariantsNewVcf.join(LohBlocks))
        // Annotate integrated SVs
        AnnotSvs = PSV_ANNOTSVS(RefWithExtras
            .combine(Metadata.map{ [it[5], it[0]] }, by: 0)
            .map{ [it[5], it[1], it[2], it[3], it[4]] }
            .join(IntegratedSvs)
        )
        // Merge LOH blocks with the GFF
        MERGE_GFFLOH(ReferenceGff
            .combine(Metadata.map{ [it[5], it[0]] }, by: 0)
            .map{ [it[2], it[0], it[1]] }
            .join(LohBlocks)
        )
}