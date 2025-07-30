include {withLog; stubLog} from '../../util/logs.nf'

def buildCmd(id, indexDir, indexPrefix, fastq1, fastq2, aligner, flags, cpus) {
    def index = "${indexDir}/${indexPrefix}"
    def align = ""
    
    if (aligner in ["bwa-mem2", "bwa", ]) {
        // Use flags.minMapq if provided as default, or override with -T if present in bwaFlags.
        def bwaFlags = flags.bwaFlags instanceof List ? flags.bwaFlags : [flags.bwaFlags]
        flags = flags + [bwaFlags: bwaFlags]

        if (flags.minMapq instanceof Integer && !flags.bwaFlags.any{it.contains("-T")}) {
            flags += [bwaFlags: flags.bwaFlags + ["-T ${flags.minMapq}"]]
        }
        def bwaFlagsArgs = flags.bwaFlags.collect{"'${it}'"}.join(" ")

        align = "${aligner} mem -t ${cpus} ${bwaFlagsArgs} '${index}' '${fastq1}' '${fastq2}'"
    } else if (aligner in ["bwameth", "bwameth-mem2"]) {
        align = "bwameth.py --threads ${task.cpus} --reference '${index}' '--do-not-penalize-chimeras' '${fastq1}' '${fastq2}'"
    } else {
        assert false, "Aligner was ${aligner}, but valid options are limited to bwa, bwa-mem2, bwameth, and bwameth-mem2."
    }
    
    def cmd = "${align} | samtools view -b -o '${id}.bam'"
    def inputMap = [id: id, fastq1: fastq1, fastq2: fastq2, aligner: aligner, index: index, flags: flags]
    def logMap = [task: "BwaAlignMates", output: "${id}.bam", input: inputMap]
    return [cmd, logMap]
}

process BwaAlignMates {
    publishDir params.general.publish.align ? params.general.publish.align : "results",
               saveAs: {params.general.publish.align ? it : null},
               mode: params.general.publish.mode
    
    label 'whenLocal_allConsuming'
    label 'align'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.alignmentContainer
    
    input:
    tuple val(id), path(indexDir), val(indexPrefix), path(fastq1), path(fastq2), val(aligner), val(flags)

    output:
    tuple val(id), path("${id}.bam")

    shell:
    (cmd, logMap) = buildCmd(id, indexDir, indexPrefix, fastq1, fastq2, aligner, flags, task.cpus)
    withLog(cmd, logMap)

    stub:
    stub = "touch '${id}.bam'"
    (cmd, logMap) = buildCmd(id, indexDir, indexPrefix, fastq1, fastq2, aligner, flags, task.cpus)
    stubLog(stub, cmd, logMap)
}