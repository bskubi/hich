include {createCompositeStrategy; filterSamplesByStrategy; columns; skip} from './extraops.nf'

process HicrepCombos{
    publishDir "results/hicrep",
               mode: params.general.publish.mode
    //container "bskubi/hich:latest"

    input:
    tuple val(planName), path(mcools), val(resolutions), val(chroms), val(exclude), val(chromFilter), val(h), val(dBPMax), val(bDownSample)

    output:
    path(outputFile)

    shell:
    outputFile = "${planName}.tsv"

    cmd = ["hich hicrep",
           resolutions ? "--resolutions ${resolutions.join(",")}" : "",
           chroms ? "--chroms ${chroms.join(",")}" : "",
           exclude ? "--exclude ${exclude.join(",")}" : "",
           chromFilter ? "--chrom-filter '${chromFilter}'" : "",
           "--h ${h.join(",")}",
           "--d-bp-max ${dBPMax.join(",")}",
           "--b-downsample ${bDownSample.join(",")}",
           "--output ${outputFile}",
           "${mcools.join(" ")}"].join(" ")
    print(cmd)
    cmd

    stub:
    outputFile = "${planName}.tsv"
    "touch ${outputFile}"
}

workflow Hicrep {
    take:
    samples

    main:

    if (!skip("hicrep")) {
        params.hicrep.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)
            filterSamplesByStrategy(samples, strategy)
                | collect
                | filter{it.size() >= 2}
                | map{columns(it, ["dropNull":true])}
                | map{
                    samples ->
                    tuple(planName,
                        samples.latestMatrix,
                        analysisPlan.resolutions,
                        analysisPlan.chroms,
                        analysisPlan.exclude,
                        analysisPlan.chromFilter,
                        analysisPlan.h,
                        analysisPlan.dBPMax,
                        analysisPlan.bDownSample)
                }
                | HicrepCombos
        }
    }


    emit:
    samples
}