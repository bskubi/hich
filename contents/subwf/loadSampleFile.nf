workflow LoadSampleFile {
    main:
    sampleFile = params.general.sampleFile
    print(sampleFile)
    samples = channel.fromPath(sampleFile.filename, checkIfExists: true)
        | splitCsv(header: true, sep: sampleFile.sep)
        

    emit:
    samples
}