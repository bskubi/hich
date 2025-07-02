include {sampleFromFastqPairs} from '../../util/cli.nf'

workflow ParseSamplesFromSRA {
    take:
    samples

    main:

    if (params.containsKey("samplesFromSRA")) {
        channel.fromSRA(params.samplesFromSRA)
            | map {
                header, files ->

                if (files.size() == 2) {
                    sample = sampleFromFastqPairs(files)
                }
                sample
            }
            | concat(samples)
            | set{samples}
    }

    emit:
    samples
}