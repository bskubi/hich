include {parsePattern} from '../../util/cli.nf'

workflow ParseFastqPairs {
    take:
    samples

    main:
    if (params.containsKey("fastqPairs")) {
        channel.fromFilePairs(params.fastqPairs, checkIfExists: true)
            | map{
                header, files ->

                sample = [fastq1: files[0], fastq2: files[1], datatype: "fastq"]

                if (params.containsKey("paramsFromPath")) {
                    fileParams = parsePattern(files[0].toString(), params.paramsFromPath)
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
