workflow LoadSampleFile {
    main:
    sampleFile = params.general.sampleFile
    samples = channel.fromPath(sampleFile.filename, checkIfExists: true)
        | splitCsv(header: true, sep: sampleFile.sep)
        | map {hmap ->
            if (hmap.get("id") == null || hmap.get("id").toString().trim().length() == 0) {
                hmap.id = "${hmap.condition}_${hmap.biorep}_${hmap.techrep}".toString()
            }
            hmap
        }

    emit:
    samples
}