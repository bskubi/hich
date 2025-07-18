include {smartFileFinder} from '../../util/files.nf'

workflow ParseSampleFile {
    take:
    samples

    main:

    if (params.containsKey("sampleFile")) {
        sampleFile = smartFileFinder(params.sampleFile).toString()

        channel.fromPath(sampleFile, checkIfExists: true)
            | splitCsv(header: true, sep: params.sampleFileSep)
            | set {newSamples}

        newSamples
            | count
            | map {count -> assert count > 0 : "Sample file is empty." }

        newSamples
            | concat(samples)
            | set{samples}
    }

    emit:
    samples
}