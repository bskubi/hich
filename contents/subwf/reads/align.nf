include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {withLog; stubLog} from '../util/logs.nf'

process BwaAlignMates {
    publishDir params.general.publish.align ? params.general.publish.align : "results",
               saveAs: {params.general.publish.align ? it : null},
               mode: params.general.publish.mode
    
    label 'whenLocal_allConsuming'
    label 'align'
    conda "$projectDir/env/dev_env.yml"
    
    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.
    
    input:
    tuple val(id), path(indexDir), val(indexPrefix), path(fastq1), path(fastq2), val(aligner), val(flags)

    output:
    tuple val(id), path("${id}.bam")

    shell:
    flags.bwaFlags = flags.bwaFlags.getClass() in [List, ArrayList] ? flags.bwaFlags : [flags.bwaFlags]

    // Use flags.minMapq if provided as default, or override with -T if present in bwaFlags.
    if (flags.minMapq instanceof Integer && !flags.bwaFlags.any{it.contains("-T")}) {
        flags.bwaFlags += ["-T ${flags.minMapq}"]
    }
    bwaFlags = flags.bwaFlags.collect{"'${it}'"}.join(" ")



    cmd = ""
    logMap = [:]
    if (aligner in ["bwa-mem2", "bwa", ]) {
        align = "${aligner} mem -t ${task.cpus} ${bwaFlags} '${indexDir}/${indexPrefix}' '${fastq1}' '${fastq2}'"
        tobam = "samtools view -b -o '${id}.bam'"
        logMap = [task: "BwaAlignMates", output: "${id}.bam", input: [id: id, fastq1: fastq1, fastq2: fastq2, aligner: aligner, index: "${indexDir}/${indexPrefix}", flags: flags]]
        cmd = "${align} | ${tobam}"
    } else if (aligner == "bsbolt") {
        logMap = [task: "BwaAlignMates", output: "${id}.bam", input: [id: id, fastq1: fastq1, fastq2: fastq2, aligner: aligner, index: indexDir, flags: flags]]
        cmd = "python3 -m bsbolt Align -t ${task.cpus} -OT ${task.cpus} -O '${id}' -DB '${indexDir}' '${bwaFlags}' -F1 '${fastq1}' -F2 '${fastq2}'"
    }
    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}.bam'"
    flags.bwaFlags = flags.bwaFlags.getClass() in [List, ArrayList] ? flags.bwaFlags : [flags.bwaFlags]

    // Use flags.minMapq if provided as default, or override with -T if present in bwaFlags.
    if (flags.minMapq instanceof Integer && !flags.bwaFlags.any{it.contains("-T")}) {
        flags.bwaFlags += ["-T ${flags.minMapq}"]
    }
    bwaFlags = flags.bwaFlags.collect{"'${it}'"}.join(" ")



    cmd = ""
    logMap = [:]
    if (aligner in ["bwa-mem2", "bwa", ]) {
        align = "${aligner} mem -t ${task.cpus} ${bwaFlags} '${indexDir}/${indexPrefix}' '${fastq1}' '${fastq2}'"
        tobam = "samtools view -b -o '${id}.bam'"
        logMap = [task: "BwaAlignMates", output: "${id}.bam", input: [id: id, fastq1: fastq1, fastq2: fastq2, aligner: aligner, index: "${indexDir}/${indexPrefix}", flags: flags]]
        cmd = "${align} | ${tobam}"
    } else if (aligner == "bsbolt") {
        logMap = [task: "BwaAlignMates", output: "${id}.bam", input: [id: id, fastq1: fastq1, fastq2: fastq2, aligner: aligner, index: indexDir, flags: flags]]
        cmd = "python3 -m bsbolt Align -t ${task.cpus} -OT ${task.cpus} -O '${id}' -DB '${indexDir}' '${bwaFlags}' -F1 '${fastq1}' -F2 '${fastq2}'"
    }

    stubLog(stub, cmd, logMap)
}

process BwaAlignSingle {
    publishDir params.general.publish.align ? params.general.publish.align : "results",
               saveAs: {params.general.publish.align ? it : null},
               mode: params.general.publish.mode
    
    label 'whenLocal_allConsuming'
    label 'align'
    
    // NOTE: Alignment speed is trivially parallelizeable and does not benefit
    // from running alignment in parallel multiple files at once. Each instance
    // of bwa-mem2 uses about 15 gigs of memory. For these two reasons we 
    // tell nextflow to run one alignment process at a time with maxForks 1.
    
    input:
    tuple val(id), path(indexDir), val(indexPrefix), path(fastq1), val(aligner), val(flags)

    output:
    tuple val(id), path("${id}.bam")

    shell:
    cmd = ""
 
    flags.bwaFlags = flags.bwaFlags.getClass() in [List, ArrayList] ? flags.bwaFlags : [flags.bwaFlags]

    // Use flags.minMapq if provided as default, or override with -T if present in bwaFlags.
    if (flags.minMapq && !flags.bwaFlags.any{it.contains("-T")}) {
        flags.bwaFlags += ["'-T ${flags.minMapq}'"]
    }
    bwaFlags = flags.bwaFlags.collect{"'${it}'"}.join(" ")

    logMap = [:]
    if (aligner in ["bwa-mem2", "bwa", ]) {
        align = "${aligner} mem -t ${task.cpus} '${bwaFlags}'' '${indexDir}/${indexPrefix}' '${fastq1}'"
        tobam = "samtools view -b -o '${id}.bam'"
        logMap = [task: "BwaAlignSingle", output: "${id}.bam", input: [id: id, fastq1: fastq1, aligner: aligner, index: "${indexDir}/${indexPrefix}", flags: flags]]
        cmd = "${align} | ${tobam}"
    } else if (aligner == "bsbolt") {
        logMap = [task: "BwaAlignSingle", output: "${id}.bam", input: [id: id, fastq1: fastq1, aligner: aligner, index: indexDir, flags: flags]]
        cmd = "python3 -m bsbolt Align -t ${task.cpus} -OT ${task.cpus} -O '${id}' -DB '${indexDir}' '${bwaFlags}' -F1 '${fastq1}'"
    }
    
    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}.bam'"
    cmd = ""
 
    flags.bwaFlags = flags.bwaFlags ?: []

    // Use flags.minMapq if provided as default, or override with -T if present in bwaFlags.
    if (flags.minMapq && !flags.bwaFlags.any{it.contains("-T")}) {
        flags.bwaFlags += ["'-T ${flags.minMapq}'"]
    }
    bwaFlags = flags.bwaFlags.collect{"'${it}'"}.join(" ")

    logMap = [:]
    if (aligner in ["bwa-mem2", "bwa", ]) {
        align = "${aligner} mem -t ${task.cpus} '${bwaFlags}'' '${indexDir}/${indexPrefix}' '${fastq1}'"
        tobam = "samtools view -b -o '${id}.bam'"
        logMap = [task: task, output: "${id}.bam", input: [id: id, fastq1: fastq1, aligner: aligner, index: "${indexDir}/${indexPrefix}", flags: flags]]
        cmd = "${align} | ${tobam}"
    } else if (aligner == "bsbolt") {
        logMap = [task: task, output: "${id}.bam", input: [id: id, fastq1: fastq1, aligner: aligner, index: indexDir, flags: flags]]
        cmd = "python3 -m bsbolt Align -t ${task.cpus} -OT ${task.cpus} -O '${id}' -DB '${indexDir}' '${bwaFlags}' -F1 '${fastq1}'"
    }
    
    stubLog(stub, cmd, logMap)
}

workflow Align {
    take:
    samples

    main:

    if (!skip("align")) {
        samples
            | filter{it.datatype == "fastq" && it.fastq1 && it.fastq2 && it.fastq1 != it.fastq2}
            | map{tuple(it.id, it.alignerIndexDir, it.alignerIndexPrefix, it.fastq1, it.fastq2, it.aligner, it.subMap("bwaFlags", "minMapq"))}
            | BwaAlignMates
            | map{[id:it[0], sambam:it[1], latest:it[1], latestSambam:it[1]]}
            | set{matesResult}
        keyUpdate(samples, matesResult, "id") | set{samples}

        samples
            | filter{it.datatype == "fastq" && it.fastq1 && (!it.fastq2 || it.fastq1 == it.fastq2)}
            | map{tuple(it.id, it.alignerIndexDir, it.alignerIndexPrefix, it.fastq1, it.aligner, it.subMap("bwaFlags", "minMapq"))}
            | BwaAlignSingle
            | map{[id:it[0], sambam:it[1], latest:it[1], latestSambam:it[1]]}
            | set{singleResult}
        keyUpdate(samples, singleResult, "id") | set{samples}
    }


    samples = emptyOnLastStep("align", samples)

    emit:
    samples
}

