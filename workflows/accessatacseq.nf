/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { INPUT_CHECK    } from '../subworkflows/local/input_check'
include { FASTQ_ALIGN_DEDUP_BWAMETH } from '../subworkflows/nf-core/fastq_align_dedup_bwameth/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { CAT_FASTQ                 } from '../modules/nf-core/cat/fastq/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
include { FASTQC                    } from '../modules/nf-core/fastqc/main'
include { TRIMGALORE                } from '../modules/nf-core/trimgalore/main'

include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_accessatacseq_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow ACCESSATACSEQ {

    take:
    samplesheet         // channel: [ path(samplesheet.csv) ]
    ch_fasta            // channel: [ path(fasta) ] // Input FASTA file for alignment
    ch_fasta_index      // channel: [ path(fasta index)     ]
    ch_bwameth_index    // channel: [ path(bwameth index)   ]
    
    main:
    // ch_fastq         = Channel.empty()
    // ch_fastqc_html   = Channel.empty()
    // ch_fastqc_zip    = Channel.empty()
    // ch_reads         = Channel.empty()
    // ch_bam           = Channel.empty()
    // ch_bai           = Channel.empty()
    // ch_qualimap      = Channel.empty()
    // ch_preseq        = Channel.empty()
    ch_multiqc_files = Channel.empty()
    ch_versions      = Channel.empty()

    
    // Branch channels from input samplesheet channel
    ch_samplesheet = samplesheet.branch { meta, fastqs ->
                            single  : fastqs.size() == 1
                                return [ meta, fastqs.flatten() ]
                            multiple: fastqs.size() > 1
                                return [ meta, fastqs.flatten() ]
                        }

    // samplesheet.map { meta, files ->
    //     println "meta.id = ${meta.id}"
    //     tuple(meta, files)
    // }

    
    // MODULE: Concatenate FastQ files from same sample if required
    CAT_FASTQ (
        ch_samplesheet.multiple
    )
    ch_fastq    = CAT_FASTQ.out.reads.mix(ch_samplesheet.single)
    ch_versions = ch_versions.mix(CAT_FASTQ.out.versions.first())

    //
    // MODULE: Run FastQC
    //
    FASTQC (
        ch_fastq
    )
    ch_fastqc_html   = FASTQC.out.html
    ch_fastqc_zip    = FASTQC.out.zip
    ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{ meta, zip -> zip })
    ch_versions      = ch_versions.mix(FASTQC.out.versions.first())

    
    // MODULE: Run TrimGalore!
    if (!params.skip_trimming) {
        TRIMGALORE(
            ch_fastq
        )
        ch_reads    = TRIMGALORE.out.reads
        ch_versions = ch_versions.mix(TRIMGALORE.out.versions.first())
    } else {
        ch_reads    = ch_fastq
    }


    //
    // SUBWORKFLOW: Alignment with BWA-METH & BAM QC
    //
    FASTQ_ALIGN_DEDUP_BWAMETH (
            ch_reads,
            ch_fasta,
            ch_fasta_index.map{ index -> [ [:], index ]},
            ch_bwameth_index,
            false, 
    )
    ch_versions = ch_versions.mix(FASTQ_ALIGN_DEDUP_BWAMETH.out.versions.unique{ it.baseName })

    // // Collect FastQC outputs from the subworkflow
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     FASTQ_ALIGN_DEDUP_BWAMETH.out.multiqc
    // )

    // // Collect outputs from the subworkflow
    // ch_bam = FASTQ_ALIGN_DEDUP_BWAMETH.out.bam
    // ch_bai = FASTQ_ALIGN_DEDUP_BWAMETH.out.bai

    //
    // Collate and save software versions
    //
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_'  +  'accessatacseq_software_'  + 'mqc_'  + 'versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }


    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    
    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList(),
        [],
        []
    )

    emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
