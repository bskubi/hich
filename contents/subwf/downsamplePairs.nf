include {transpack; isTechrep; isBiorep; isCondition; emptyOnLastStep} from './extraops.nf'

/*
    1. Compute stats over the pairs file with hich stats as --orig-stats
        This is a single process call and attaches the stats file to the sample

    2. Set strata and aggregate the stats to obtain target stats as --target-stats
        This is possibly multiple process calls (?) and attaches the stats file to the sample

    3. Use hich downsample with --orig_stats and --target_stats
        This is a single function call with transpack that specifies the new id, similar to merge
    
    There was an idea to modify transpack so that you could pass in a groupTuple output
    and have it collect the parameters from each sample in each group and pass it in. This
    could be used for the merge steps as well as here.
        We pass in one or more channels to transpack/transact
        Transpack would just pass through the first channel as normal, no modifications necessary.
        Transact would be modified during the map step.
            If the channel items are HashMaps, do things as normal
            If the channel items are two-element tuples [groupBy, [HashMaps]]
            then:
                Collect the values from each HashMap associated with the input arguments
                Pass to the function
*/

process HichDownsample {
    input:
    tuple val(id), path(pairs), path(orig_stats), path(target_stats)

    output:
    tuple val(id), path("${id}_downsampled.pairs.gz")

    shell:
    "hich downsample --orig_stats ${orig_stats} --target_stats ${target_stats} ${pairs} ${id}_downsampled.pairs.gz"

    stub:
    "touch ${id}_downsampled.pairs.gz"
}

process HichStats {
    input:
    tuple val(id), path(pairs), val(hich_stats_params)

    output:
    tuple val(id), path("${id}.stats.tsv")

    shell:
    hich_stats_params = hich_stats_params.join(" ")
    "hich stats ${hich_stats_params} ${pairs} ${id}.stats.tsv"
}

process HichStatsAggregateSubgroups {
    input:
    tuple val(subgroup), val(ids), path(stats), val(statsAggregateParams), val(stats_outputs)

    output:
    tuple val(subgroup), val(ids), path(stats_outputs)

    shell:
    def toAggregate = (pairs.getClass() == nextflow.processor.TaskPath) ? samples : samples.join(" ")
    "hich stats-aggregate ${statsAggregateParams} ${toAggregate}"
}

process HichStatsAggregateAll {
    input:
    tuple val(ids), path(stats), val(statsAggregateParams), val(stats_outputs)

    output:
    tuple val(ids), path(stats_outputs)

    shell:
    def toAggregate = (pairs.getClass() == nextflow.processor.TaskPath) ? samples : samples.join(" ")
    "hich stats-aggregate ${statsAggregateParams} ${toAggregate}"
}


workflow Downsample {
    take:
    samples
    sampleType

    main:
    // Note: we need a way of skipping downsampling where it's not desired.
    // Get the sample traits to group by ("groups") and the sample trait name
    // to save the resulting pairs file as ("downsample_techreps_pairs")
    switch(sampleType) {
        case "DownsampleTechreps":
            groups = ["condition", "biorep"]
            downsampled_pairs = "downsample_techreps_pairs"
            orig_stats = "techreps_stats"
            target_stats = "downsample_techreps_target"
            break
        case "DownsampleBioreps":
            groups = ["condition"]
            downsampled_pairs = "downsample_bioreps_pairs"
            orig_stats = "bioreps_stats"
            target_stats = "downsample_bioreps_target"
            break
        case "DownsampleConditions":
            groups = []
            downsampled_pairs = "downsample_conditions_pairs"
            orig_stats = "conditions_stats"
            target_stats = "downsample_techreps_target"
            break
        default:
            error "Invalid merge sampleType '${sampleType}'"
    }

    // First, filter the samples for only techreps, bioreps, or conditions as appropriate.
    samples
        | filter{
            switch(sampleType) {
                case "DownsampleTechreps":
                    isTechrep(it)
                    break
                case "DownsampleBioreps":
                    isBiorep(it)
                    break
                case "DownsampleConditions":
                    isCondition(it)
                    break
                default:
                    error "Invalid merge sampleType '${sampleType}'"
            }
        }
        | map{it[orig_stats]}
        | set{to_downsample}

    // Second, call stats on to_downsample with transpack
    // to_downsample = transpack(
    //     HichStats,
    //     [to_downsample],
    //     ["id", "latest", "hichStatsParams"],
    //     ["id", orig_stats],
    //     [],
    //     "id"
    // )

    // Then put the traits ("condition", "stats", "none") as the first item in a tuple with the items themselves as the second
    // Then do a groupTuple on the first item to group the filtered samples with their similars
    // (other techreps in the same condition + biorep). This will have the form [[grouping], [list of HashMap items]]
    to_downsample
        | map{tuple(it.subMap(*groups), it)}
        | groupTuple
    
    
    // Call stats-aggregate and store as target distribution
    // This is where I'm still a bit confused about implementation.
    // At this point, we've filtered for the appropriate condition and split to subgroups.
    // However, as currently written, stats-aggregate only offers the options:
    // --to_group_mean
    // --to_group_min
    // --to_count
    // If we pass the individual groups to the process that gives us the ability to do subgroup
    // filtering and obtain output stats files, but it doesn't let us do additional processing
    // on the entire group. So I think what we'd need to do here is something like:
    // call stats-aggregate once per subgroup using a given prefix
    // construct the step 1 output stats file names based on the prefix
    // then call stats-aggregate again on the entire group using the step 1
    // construct the step 2 output stats file names based ont the prefix
    // then return the step 2 output stats file names as the targets
    //
    // The problem here is that to stage the files appropriately we have to pass them in as paths
    // but we can't do that with a list of lists.
    // We could do that by nesting the ids and passing in all the paths
    // then we'd construct the pipeline in the process
    // Potentially we could do the whole thing in a single process.


    // Call downsample and store as latest + "downsampled"
    // samples = transpack(
    //     HichStats,
    //     [to_downsample, samples],
    //     ["id", "latest", orig_stats, target_stats],
    //     ["id", downsampled_pairs],
    //     ["latest":downsampled_pairs],
    //     "id"
    // )

    emit:
    samples
}

workflow DownsampleConditions {
    take:
    samples

    main:
    samples = Downsample(samples, "DownsampleConditions")

    samples = emptyOnLastStep("DownsampleConditions", samples)

    emit:
    samples
}


workflow DownsampleBioreps {
    take:
    samples

    main:
    samples = Downsample(samples, "DownsampleBioreps")

    samples = emptyOnLastStep("DownsampleBioreps", samples)


    emit:
    samples
}

workflow DownsampleTechreps {
    take:
    samples

    main:
    samples = Downsample(samples, "DownsampleTechreps")

    samples = emptyOnLastStep("DownsampleTechreps", samples)
    
    emit:
    samples
}