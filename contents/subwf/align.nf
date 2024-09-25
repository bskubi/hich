include {transpack; emptyOnLastStep; updateChannel} from './extraops.nf'

process BwaAlign {
    publishDir params.general.publish.align ? params.general.publish.align : "results",
               saveAs: {params.general.publish.align ? it : null},
               mode: params.general.publish.mode

    conda "bioconda::bwa bioconda::bwa-mem2 bioconda::samtools"
    container "bskubi/hich:latest"
    
    label 'whenLocal_allConsuming'
    label 'align'
    
    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.

    // NOTE 2: I need to decide on how to work with the possibility that users
    // will run multiple aligners (bwa-mem2, bwa, bsbolt, etc). We should still
    // only run one at a time.

    input:
    tuple val(id), path(index_dir), val(index_prefix), path(fastq1), path(fastq2), val(aligner), val(bwaFlags)

    output:
    tuple val(id), path("${id}.bam")

    shell:
    align = ""
    if (aligner in ["bwa-mem2", "bwa"]) {
        align = "${aligner} mem -t ${task.cpus} ${bwaFlags} ${index_dir}/${index_prefix} ${fastq1} ${fastq2}"
    }
    
    tobam = "samtools view -b -o ${id}.bam"

    "${align} | ${tobam}"

    stub:
    "touch ${id}.bam"
}

workflow Align {
    take:
    samples

    main:

    samples | branch{yes: it.datatype == "fastq"; no: true} | set {run}
    run.yes
        | map{tuple(it.id, it.alignerIndexDir, it.alignerIndexPrefix, it.fastq1, it.fastq2, it.aligner, it.bwaFlags)}
        | BwaAlign
        | map{[id:it[0], sambam:it[1], latest:it[1], latestSambam:it[1]]}
        | set{runResult}
    updateChannel(run.yes, runResult) | concat(run.no) | set{samples}

    samples = emptyOnLastStep("Align", samples)

    emit:
    samples
}