

workflow Split {
    take:
    levelSamples
    levelParams

    main:
    // levelParams gets converted to a DataflowVariable when emitted from the previous workflow,
    // so extract the LinkedHashMap value so we can use it more conveniently
    levelParams = levelParams.value

    levelSamples | branch {
        YES: (
            // Select samples that have an aggregateProfileName and that include
            // a value of either downsampleToMeanDistribution or downsampletoSize
            // that requires processing.
            it.aggregateProfileName != null 
            && (
                it.get(levelParams.downsampleToMeanDistribution) != null
                || !(it.get(levelParams.downsampleToSize) in [1.0, null])
            )
        )
        NO: true
    } | set{downsample}

    // Compute the number of reads downsampled from
    downsample.YES
        | map{tuple(it.id, it.latestPairs, it.get(levelParams.readConjuncts), it.get(levelParams.cisStrata))}
        | HichStats
        | map{[id:it[0], (levelParams.downsampleStatsFrom):it[1]]}
        | set{downsampleStatsFromResult}

    // Join the results to the sample map
    pack(downsample.YES, downsampleStatsFromResult) | set {fromStatsCalculated}

    // Compute the number of reads downsampled to
    // These are converted to columnar format as input to HichStatsAggregate,
    // then converted back to rows format and flattened into channel elements
    groupHashMap(fromStatsCalculated, [levelParams.downsampleToMeanDistribution, 'aggregateProfileName'])
        | map{columns(it, ["dropNull":true])}
        | map{tuple(it.id, it.get(levelParams.downsampleStatsFrom), it.outlier)}
        | HichStatsAggregate
        | map{[id:it[0], (levelParams.downsampleStatsTo):it[1]]}
        | map{rows(it)}
        | flatten
        | set {downsampleToGroupMinResult}

    // Join the results to the sample map
    pack(fromStatsCalculated, downsampleToGroupMinResult) | set{toGroupMinResult}

    // Downsample from the number of reads in each bin in levelParams.downsampleStatsFrom to levelParams.downsampleStatsTo
    toGroupMinResult
        | map{
            tuple(it.id, it.latestPairs, it.get(levelParams.downsampleStatsFrom), it.get(levelParams.downsampleStatsTo),
              it.get(levelParams.readConjuncts), it.get(levelParams.cisStrata), it.get(levelParams.downsampleToSize))}
        | HichDownsample
        | map{[id:it[0], (levelParams.downsamplePairs):it[1], latest:it[1], latestPairs:it[1]]}
        | set{downsampledResult}

    // Join the results to the sample map
    pack(toGroupMinResult, downsampledResult) | set {downsampled}

    // Add back in the non-downsampled samples
    downsampled | concat(downsample.NO) | set{result}

    emit:
    result
    levelParams
}
