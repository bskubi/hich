include {transpack; isTechrep; isBiorep; isCondition; emptyOnLastStep; transposeHashMapList} from './extraops.nf'

/*
    1. Compute stats over the pairs file with hich stats as --orig-stats
        This is a single process call and attaches the stats file to the sample

    2. Set strata and aggregate the stats to obtain target stats as --target-stats
        This is possibly multiple process calls (?) and attaches the stats file to the sample

    3. Use hich downsample with --orig_stats and --target_stats
        This is a single function call with transpack that specifies the new id, similar to merge
    
    We can now use transposeHashMapList to get the input stats file for a given group in column order
    and use the options = [transpose:true] flag to get multiple outputs.

    WHY have separate workflows for downsampling techreps, bioreps, and conditions?
        BECAUSE if we want to let the user provide equal techrep input sizes when merging to bioreps,
        or equal biorep input sizes when merging to conditions, we need to have downsampled these
        groups before the merge happens.

    FEATURES
        Let the user set samples as outliers and declare whether that carriers over during the merge
    
    The basic feature set we've been trying to build is something like:
        Let the user use outlier-aware stratified downsampling to the minimum group size
        Let the user downsample to a specific integer or % size

        One way to approach this would be to give options for:

        
        Define 1+ sets of conjuncts -> match all techreps within a biorep (yes, no, both) -> then downsample to a specific size
        Define 1+ sets of conjuncts -> match all bioreps within a condition (yes, no both) -> then downsample to a specific size
        Define 1+ sets of conjuncts -> match all conditions (yes no both) -> then downsample to a specific size

        Of course the user might want all combinations of these, or only a subset.
        Usually we just deliver all combinations. However this could lead to very long runtimes.

        One option is to deliver all combinations, and require the user to define small subsets of the combinations
        and rerun if they don't want all possible combinations.

        Another is to require the user to specify each combination individually.
        This is something the user can also do with Hich directly. I think if they have specific requirements they can just do that.

Example with defining each combination individually:

coverageControl {
    all {}

    downsampleByHalf {
        Any samples matching any of the HashMaps in outliers will be treated as such
        outliers = [[condition: Mock, biorep: 2, techrep: [1, 2, 3]]]

        For all techreps in the same condition + biorep, downsample and label with "downsampleByHalf"
        as their coverageControl profile. Samples are merged based on their coverageControl profile.
        The default coverageControl profile is 'all', which means no downsampling.

        Plausibly users might want to not process the 'all', perhaps not merging it and not turning it into a contact matrix.
        This suggests an overall feature of skipping specific steps for specific samples. We were already going to implement
        a step skipping feature. We could make both the lastStep and the skipStep features sample specific. This would let the user
        use --lastStep CoverageControl [coverageProfile: all].

        One tricky aspect is figuring out how to pass these sorts of arguments to Nextflow. Possibly this would just have to be done in the config.

        The other possible issue here though is that you could have a coverage profile that, let's say, doesn't downsample techreps or bioreps,
        just conditions. That means 'all' techreps need to be merged. But if you set --lastStep CoverageControl [coverageProfile: all] you wouldn't have
        the 'all' bioreps.

        There's an inefficiency here which is that if two downsampling profiles have the same starting steps, they have to be done repeatedly (I think?)
        unless we can set up the process so that the inputs would be identical in either case so it just uses the cached results.

        Another possibility is that we bundle the downsampling AND merge steps into one merge profile.
        That way we can define the needed inputs for both steps.

        Here, we have no downsampling parameters over the techreps so we just merge the 'all' samples to bioreps.
        mergeTechrepsToBioreps: true
        
        Because we define ways of downsampling the bioreps, they do get downsampled after being merged. However, since there is no
        mergeBiorepsToConditions, they don't get merged.
        biorepsConjuncts = ['record.chr1 record.chr2 stratum']
        matchCoverage: "subgroup" or "all"
        downsampleBioreps: .5
    }
}

I think this looks very good as an interface. Explicit, controlled, and sparse, which I think is the default mode most people will want.

So then, we can combine the downsample and merge steps.
We only need to produce stats IF we are downsampling.


Downsampling: true -> Collect stats using *Conjuncts
    matchCoverage: "subgroup" -> aggregate subgroups  "min"     We can simply have this be a process
                   "all"      -> aggregate all groups "min"     And this be a process
    downsample*: true         -> aggregate .5                   And this be a process
    downsample                                                  And this be a process
Merge: true -> Merge                                            And this be a process

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