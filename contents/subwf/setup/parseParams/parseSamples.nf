include {parsePattern} from '../../util/cli.nf'
include {formatFromExtension} from '../../util/files.nf'

workflow ParseSamples {
    take:
    samples

    main:

    if (params.containsKey("samples")) {
        channel.fromPath(params.samples)
            | map {
                file ->
                format = formatFromExtension(file.toString())
                sample = [(format):file]

                switch (format) {
                    case "sambam":
                        sample += [latestSambam: file, datatype: "sambam"]
                        break
                    case "pairs":
                        sample += [latestPairs: file, datatype: "pairs"]
                        break
                    case "mcool":
                        sample += [latestMatrix: file, datatype: "matrix"]
                        break
                    case "hic":
                        sample += [latestMatrix: file, datatype: "matrix"]
                        break
                }



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