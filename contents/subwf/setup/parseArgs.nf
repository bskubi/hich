include {withLog; stubLog; emptyOnLastStep; parsePattern; datatypeFromExtension} from '../extraops.nf'

def sampleFromFastqPairs(files) {
    file1 = files[0]
    file2 = files[1]

    bothFastq = datatypeFromExtension(file1) == "fastq" && datatypeFromExtension(file2) == "fastq"
    assert bothFastq, "--fastqPairs grouped ${file1} and ${file2}, but these do not contain a .fastq or .fq extension as expected."
    sample = [fastq1:file1, fastq2: file2, datatype: "fastq"]

    if (params.containsKey("paramsFromPath")) {
        f1Params = parsePattern(file1.toString(), params.paramsFromPath)
        f2Params = parsePattern(file2.toString(), params.paramsFromPath)
        sameParams = f1Params == f2Params
        assert sameParams, "--paramsFromPath yielded different params for ${file1} and ${file2}:\n${f1Params}\n${f2Params}"
        sample += f1Params
    }
    else if (params.containsKey("paramsFromFilename")) {
        f1Params = parsePattern(file1.name.toString(), params.paramsFromFilename)
        f2Params = parsePattern(file2.name.toString(), params.paramsFromFilename)
        sameParams = f1Params == f2Params
        assert sameParams, "--paramsFromFilename yielded different params for ${file1} and ${file2}:\n${f1Params}\n${f2Params}"
        sample += f1Params
    }
    return sample
}

workflow ParseArgs {
    main:
    samples = channel.empty()

    if (params.containsKey("sampleFile")) {
        channel.fromPath(params.sampleFile, checkIfExists: true)
            | splitCsv(header: true, sep: params.sampleFileSep)
            | set{fromSampleFile}
        samples | concat(fromSampleFile) | set{samples}
    }

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
            | set {fromFastqInterleavedChannel}
        samples | concat(fromFastqInterleavedChannel) | set{samples}
    }

    if (params.containsKey("fastqPairs")) {
        channel.fromFilePairs(params.fastqPairs)
            | map{
                header, files ->

                sample = [fastq1: files[0], fastq2: files[1], datatype: "fastq"]

                if (params.containsKey("paramsFromPath")) {
                    fileParams = parsePattern(files[0].toString(), params.paramsFromPath)
                    sample += fileParams
                }


                sample
            }
            | set{fromFilePairsChannel}
        samples | concat(fromFilePairsChannel) | set{samples}
    }

    if (params.containsKey("samples")) {
        channel.fromPath(params.samples)
            | map {
                file ->
                datatype = datatypeFromExtension(file.toString())
                sample = [(datatype):file, datatype: datatype]

                if (params.containsKey("paramsFromPath")) {
                    fileParams = parsePattern(file.toString(), params.paramsFromPath)
                    sample += fileParams
                }

                sample
            }
            | set{fromFileChannel}
        samples | concat(fromFileChannel) | set{samples}
    }

    if (params.containsKey("samplesFromSRA")) {
        channel.fromSRA(params.samplesFromSRA)
            | map {
                header, files ->

                if (files.size() == 2) {
                    sample = sampleFromFastqPairs(files)
                }
                sample
            }
            | set{fromSRAChannel}
        samples | concat(fromSRAChannel) | set{samples}
    }

    samples = emptyOnLastStep("parseArgs", samples)

    emit:
    samples
}