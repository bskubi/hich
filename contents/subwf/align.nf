include {emptyOnLastStep; pack; skip} from './extraops.nf'

process BwaAlignMates {
    publishDir params.general.publish.align ? params.general.publish.align : "results",
               saveAs: {params.general.publish.align ? it : null},
               mode: params.general.publish.mode

    conda "bioconda::bwa bioconda::bwa-mem2 bioconda::samtools"
    container params.general.hichContainer
    
    label 'whenLocal_allConsuming'
    label 'align'
    
    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.
    
    input:
    tuple val(id), path(indexDir), val(indexPrefix), path(fastq1), path(fastq2), val(aligner), val(bwaFlags)

    output:
    tuple val(id), path("${id}.bam")

    shell:
    cmd = ""
    if (aligner in ["bwa-mem2", "bwa", ]) {
        align = "${aligner} mem -t ${task.cpus} '${indexDir}/${indexPrefix}' '${bwaFlags}' '${fastq1}' '${fastq2}'"
        tobam = "samtools view -b -o ${id}.bam"
        cmd = "${align} | ${tobam}"
    } else if (aligner == "bsbolt") {
        cmd = "python3 -m bsbolt Align -t ${task.cpus} -OT ${task.cpus} -O '${id}' -DB '${indexDir}' '${bwaFlags}' -F1 '${fastq1}' -F2 '${fastq2}'"
    }
    
    cmd

    stub:
    "touch '${id}.bam'"
}

process BwaAlignSingle {
    publishDir params.general.publish.align ? params.general.publish.align : "results",
               saveAs: {params.general.publish.align ? it : null},
               mode: params.general.publish.mode

    conda "bioconda::bwa bioconda::bwa-mem2 bioconda::samtools"
    container params.general.hichContainer
    
    label 'whenLocal_allConsuming'
    label 'align'
    
    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.
    
    input:
    tuple val(id), path(indexDir), val(indexPrefix), path(fastq1), val(aligner), val(bwaFlags)

    output:
    tuple val(id), path("${id}.bam")

    shell:
    cmd = ""
    if (aligner in ["bwa-mem2", "bwa", ]) {
        align = "${aligner} mem -t ${task.cpus} '${bwaFlags}'' '${indexDir}/${indexPrefix}' '${fastq1}'"
        tobam = "samtools view -b -o ${id}.bam"
        cmd = "${align} | ${tobam}"
    } else if (aligner == "bsbolt") {
        cmd = "python3 -m bsbolt Align -t ${task.cpus} -OT ${task.cpus} -O '${id}' -DB '${indexDir}' '${bwaFlags}' -F1 '${fastq1}'"
    }
    
    cmd

    stub:
    "touch '${id}.bam'"
}

workflow Align {
    take:
    samples

    main:

    samples
        | filter{!skip("align") && it.datatype == "fastq" && it.fastq1 && it.fastq2 && it.fastq1 != it.fastq2}
        | map{tuple(it.id, it.alignerIndexDir, it.alignerIndexPrefix, it.fastq1, it.fastq2, it.aligner, it.bwaFlags)}
        | BwaAlignMates
        | map{[id:it[0], sambam:it[1], latest:it[1], latestSambam:it[1]]}
        | set{matesResult}
    pack(samples, matesResult) | set{samples}

    samples
        | filter{!skip("align") && it.datatype == "fastq" && it.fastq1 && (!it.fastq2 || it.fastq1 == it.fastq2)}
        | map{tuple(it.id, it.alignerIndexDir, it.alignerIndexPrefix, it.fastq1, it.aligner, it.bwaFlags)}
        | BwaAlignSingle
        | map{[id:it[0], sambam:it[1], latest:it[1], latestSambam:it[1]]}
        | set{singleResult}
    pack(samples, singleResult) | set{samples}

    samples = emptyOnLastStep("align", samples)

    emit:
    samples
}