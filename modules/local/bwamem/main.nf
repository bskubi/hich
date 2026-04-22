process BWAMEM {
    publishDir = "results/bwamem/"
    tag "$ID"

    input:
    tuple val(ID), path(fastq), path(indexDir), val(indexPfx), val(is3DGenome), val(isInterleaved)

    output:
    tuple val(ID), path(bam)

    script:
    // Output path
    bam = "${ID}.bam"

    flags = [
        // Number of threads
        "-t ${task.cpus}",

        // Paired-end, interleaved in single fastq file
        *(isInterleaved ? ["-p"] : []),
        
        // Flags for chimeric (Hi-C) alignment
        *(is3DGenome ? ["-SP5M"] : [])
    ].join(" ")

    // convert to space-delimited list of fastq paths
    fastq = [fastq].flatten().collect { it.toString() }.join(" ")

    // format indexDir to combine with indexPfx
    index = indexDir.toString().stripEnd("/")

    """
    bwa mem ${flags} ${index}/${indexPfx} ${fastq} | samtools view -b > ${bam}
    """

    stub:
    bam = "${ID}.bam"

    """
    touch ${bam}
    """
}