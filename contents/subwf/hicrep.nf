include {createCompositeStrategy; filterSamplesByStrategy; columns; skip; formatArg} from './extraops.nf'

process HicrepCombos{
    publishDir "results/hicrep", mode: params.general.publish.mode
    //container "bskubi/hich:latest"
    
    input:
    tuple val(planName), path(mcools), val(resolutions), val(chroms), val(exclude), val(chromFilter), val(h), val(dBPMax), val(bDownSample)

    output:
    path(outputFile)

    shell:
    outputFile = "${planName}.tsv"

    cmd = ["hich hicrep",
           formatArg("--resolutions %s", resolutions, ','),
           formatArg("--chroms %s", chroms, ','),
           formatArg("--exclude %s", exclude, ','),
           formatArg("--chrom-filter '%s'", chromFilter, ''),
           formatArg("--h %s", h, ','),
           formatArg("--d-bp-max %s", dBPMax, ','),
           formatArg("--b-downsample %s", bDownSample, ','),
           "--output ${outputFile}",
           formatArg("%s", mcools, ' ')].findAll{it}.join(" ")
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
        hicrepParameterizations = channel.empty()

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
                | concat(hicrepParameterizations)
                | set{hicrepParameterizations}
        }

        hicrepParameterizations | HicrepCombos
    }


    emit:
    samples
}