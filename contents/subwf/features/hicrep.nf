include {createCompositeStrategy; filterSamplesByStrategy; groupSamplesByStrategy} from '../util/analysisPlans.nf'
include {skip; formatArg} from '../util/cli.nf'
include {columns} from '../util/reshape.nf'
include {withLog; stubLog} from '../util/logs.nf'

process HicrepCombos{
    publishDir "results/hicrep", mode: params.general.publish.mode
    tag "$planName"

    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer
    
    input:
    tuple val(planName), path(mcools), val(resolutions), val(chroms), val(exclude), val(chromFilter), val(h), val(dBPMax), val(bDownSample)

    output:
    path(hicrep)

    shell:
    hicrep = "${planName}.tsv"
    resolution = "-r " + resolutions.join(" -r")

    cmd = ["hich matrix hicrep",
           resolution,
           formatArg("--chroms %s", chroms, ','),
           formatArg("--exclude %s", exclude, ','),
           formatArg("--chrom-filter '%s'", chromFilter, ''),
           formatArg("--h %s", h, ','),
           formatArg("--d-bp-max %s", dBPMax, ','),
           formatArg("--b-downsample %s", bDownSample, ','),
           "'${hicrep}'",
           formatArg("%s", mcools, ' ')].findAll{it}.join(" ")
    logMap = [task: "HicrepCombos", input: [planName: planName, mcools: mcools, resolutions: resolutions, chroms: chroms, exclude: exclude, chromFilter: chromFilter, dBPMax: dBPMax, bDownSample: bDownSample], 
    output: [hicrep: hicrep]]
    withLog(cmd, logMap)

    stub:
    hicrep = "${planName}.tsv"
    stub = "touch ${hicrep}"
    resolution = "-r " + resolutions.join(" -r")
    cmd = ["hich matrix hicrep",
           resolution,
           formatArg("--chroms %s", chroms, ','),
           formatArg("--exclude %s", exclude, ','),
           formatArg("--chrom-filter '%s'", chromFilter, ''),
           formatArg("--h %s", h, ','),
           formatArg("--d-bp-max %s", dBPMax, ','),
           formatArg("--b-downsample %s", bDownSample, ','),
           "'${hicrep}'",
           formatArg("%s", mcools, ' ')].findAll{it}.join(" ")
    logMap = [task: "HicrepCombos", input: [planName: planName, mcools: mcools, resolutions: resolutions, chroms: chroms, exclude: exclude, chromFilter: chromFilter, dBPMax: dBPMax, bDownSample: bDownSample], 
    output: [hicrep: hicrep]]
    stubLog(stub, cmd, logMap)

}

workflow Hicrep {
    take:
    samples

    main:

    if (!skip("hicrep")) {
        hicrepParameterizations = channel.empty()

        params.hicrep.each {
            planName, analysisPlan ->

            // Create composite strategies based on keys within params.hicrep as the sampleSelectionStrategy
            // and params.sampleSelectionStrategies as the individual sample selection strategies to draw from
            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            // Filter out samples that 
            filtered = filterSamplesByStrategy(samples, strategy)
            grouped = groupSamplesByStrategy(filtered, strategy)
            grouped
                | filter{it.size() >= 2}
                | map{columns(it, ["dropAllNull":true])}
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