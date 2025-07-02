include {parsePattern} from '../../util/cli.nf'

workflow ParseFastqInterleaved {
    take:
    samples

    main:

    if (params.containsKey("fastqInterleaved")) {
        channel.fromPath(params.fastqInterleaved, checkIfExists: true)
            | map {
                fastq ->
                sample = [fastq1: fastq, fastq2: fastq, datatype: "fastq"]

                if (params.containsKey("paramsFromPath")) {
                    fileParams = parsePattern(fastq.toString(), params.paramsFromPath)
                    sample += fileParams
                }
                sample
            }
            | concat(samples)
            | set{samples}
    }

    emit:
    samples
}