include {parsePattern} from '../../util/cli.nf'
include {datatypeFromExtension} from '../../util/files.nf'

workflow ParseSamples {
    take:
    samples

    main:

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
            | concat(samples)
            | set{samples}
    }

    emit:
    samples
}