include {updateChannel; coalesce; rows; columns; groupHashMap; isTechrep; isBiorep; isCondition; constructIdentifier} from './extraops.nf'

process HichDownsampleStatsFrom {
    input:
    tuple val(id), val(conjuncts), val(cisStrata)

    output:
    tuple val(id), path(stats)

    stub:
    stats = "${id}.from.stats.tsv"
    "touch ${stats}"
}

process HichDownsampleStatsTo {
    input:
    tuple val(ids), path(stats), val(toSizes), val(outliers)

    output:
    tuple val(ids), path(targetStats)

    shell:
    targetStats = ids.collect{"${it}.to.stats.tsv"}
    "touch ${targetStats.join(" ")}"
}

process HichDownsamplePairs {
    input:
    tuple val(id), val(conjuncts), val(cisStrata), path(statsFrom), path(statsTo)

    output:
    tuple val(id), path(downsampledPairs)

    shell:
    downsampledPairs = "${id}.downsampled.pairs.tsv"
    "touch ${downsampledPairs}"
}

process PairtoolsMerge {
    input:
    tuple val(id), path(pairs)

    output:
    tuple val(id), path(merged)

    shell:
    merged = "${id}.merged.pairs.gz"
    "touch ${merged}"
}

process PairtoolsDedup {
    input:
    tuple val(id), path(pairs), val(singleCell), val(maxMismatch), val(method)

    output:
    tuple val(id), path(deduplicated)

    shell:
    deduplicated = "${id}.dedup.pairs.gz"
    "touch ${deduplicated}"
}

/*
    Todo: Refactor to keep as DRY as possible
    Todo: Create AggregateBioreps and AggregateConditions
    Todo: Redo Hich stats, stats-aggregate, and downsample interface
    Todo: Implement stats, stats-aggregate, downsample, merge and dedup

    If you have a large block of code that's only called once, should you
    split the parts out into functions? Or just label their overall purpose
    with a comment?
*/

workflow CreateAggregatePairProfiles {
    take:
    levelSamples
    levelParams

    main:
    // Split raw samples from the current level vs. those already associated with an aggregation profile
    levelSamples | branch{
        YES: it.aggregateProfileName == null
        NO: true
    } | set{raw}

    // From the raw samples, create new samples associated with each aggregation profile
    // This may need to be redone putting the aggregateProfiles into a channel of their own and using combinations.

    
    aggregateProfiles = channel.of(params.comparisonSets.aggregateProfiles)
    | map{[profileName: it.keySet().toList(), profileParams: it.values().toList()]}
    | map{rows(it)}
    | flatten
    | combine(raw.YES)
    | map {
        profile = it[0]
        sample = it[1]
        sample += [aggregateProfileName: profile.profileName] + profile.profileParams
        sample
    }
    | map{
        downsampleToMeanDistributionGroup = it.subMap((levelParams.downsampleToMeanDistribution)) ?: null
        it + [(levelParams.downsampleToMeanDistribution): downsampleToMeanDistributionGroup] + [id: "${it.id}_${it.aggregateProfileName}"]
    }
    | concat(raw.NO)
    | set{result}


    emit:
    result
    levelParams
}

workflow DownsamplePairs {
    take:
    levelSamples
    levelParams

    main:
    // levelParams gets converted to a DataflowVariable when emitted from the previous workflow,
    // so extract the LinkedHashMap value so we can use it more conveniently
    levelParams = levelParams.value

    levelSamples | branch {
        YES: it.aggregateProfileName != null && (it.get(levelParams.downsampleToMeanDistribution) != null || it.get(levelParams.downsampleToSize) != 1.0)
        NO: true
    } | set{downsample}

    // Compute from stats
    downsample.YES
        | map{tuple(it.id, it.get(levelParams.readConjuncts), it.get(levelParams.cisStrata))}
        | HichDownsampleStatsFrom
        | map{[id:it[0], (levelParams.downsampleStatsFrom):it[1]]}
        | set{downsampleStatsFromResult}
    updateChannel(downsample.YES, downsampleStatsFromResult) | set {fromStatsCalculated}

    // Compute to stats
    groupHashMap(fromStatsCalculated, [levelParams.downsampleToMeanDistribution, 'aggregateProfileName'])
        | map{columns(it, ["dropNull":true])}
        | map{tuple(it.id, it.get(levelParams.downsampleStatsFrom), it.get(levelParams.downsampleToSize), it.outlier)}
        | HichDownsampleStatsTo
        | map{[id:it[0], (levelParams.downsampleStatsTo):it[1]]}
        | map{rows(it)}
        | flatten
        | set {downsampleStatsToResult}
    updateChannel(fromStatsCalculated, downsampleStatsToResult) | set{toStatsCalculated}
    toStatsCalculated
        | map{tuple(it.id, it.get(levelParams.readConjuncts), it.get(levelParams.cisStrata), it.get(levelParams.downsampleStatsFrom), it.get(levelParams.downsampleStatsTo))}
        | HichDownsamplePairs
        | map{[id:it[0], (levelParams.downsamplePairs):it[1], latest:it[1], latestPairs:it[1]]}
        | set{downsampledResult}
    updateChannel(toStatsCalculated, downsampledResult) | set {downsampled}

    downsampled | concat(downsample.NO) | set{result}

    emit:
    result
    levelParams
}

workflow MergePairs {
    take:
    levelSamples
    levelParams

    main:

    // levelParams gets converted to a DataflowVariable when emitted from the previous workflow,
    // so extract the LinkedHashMap value so we can use it more conveniently
    levelParams = levelParams.value

    levelSamples | branch {
        YES: it.includeInMerge != false && it.get(levelParams.doMerge)
        NO: true}
    | set{merge}

    // Merge
    groupHashMap(merge.YES, levelParams.mergeGroupIdentifiers)
        | map{coalesce(columns(it, ["dropNull":true]))}
        | map{tuple(constructIdentifier(coalesce(it, "_drop")), it.latestPairs)}
        | PairtoolsMerge
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{fromMerge}

    groupHashMap(merge.YES, levelParams.mergeGroupIdentifiers)
        | map{coalesce(columns(it, ["dropNull":true]), '_drop')}
        | map{it += [id:constructIdentifier(it)]}
        | set{inheritedMergeAttributes}
    updateChannel(inheritedMergeAttributes, fromMerge) | set{merged}

    merged | concat(levelSamples) | set{result}

    /*
    !The merge outputs need to be passed through the setup protocol
    */

    emit:
    result
    levelParams
}

workflow DeduplicatePairs {
    take:
    levelSamples
    levelParams

    main:
    
    // levelParams gets converted to a DataflowVariable when emitted from the previous workflow,
    // so extract the LinkedHashMap value so we can use it more conveniently
    levelParams = levelParams.value

    levelSamples
        | branch {
            YES: it.deduplicate && !it.alreadyDeduplicated && levelParams.levelFilter(it) && it.get(levelParams.doDedup)
            NO: true
    } | set{deduplicate}

    deduplicate.YES
        | map{tuple(it.id, it.latestPairs, it.dedupSingleCell, it.dedupMaxMismatch, it.dedupMethod)}
        | PairtoolsDedup
        | map{[id:it[0], dedupPairs:it[1], latest:it[1], latestPairs:it[1], alreadyDeduplicated:true]}
        | set{dedupResult}
    updateChannel(deduplicate.YES, dedupResult) | concat(deduplicate.NO) | set{result}

    emit:
    result
}

workflow AggregateTechreps {
    take:
    samples

    main:

    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'max',
        cisStrata: 'techrepCisStrata',
        readConjuncts: 'techrepReadConjuncts',
        downsampleToMeanDistribution: 'techrepDownsampleToMeanDistribution',
        downsampleToSize: 'techrepDownsampleToSize',
        downsampleStatsFrom: 'techrepDownsampleStatsFrom',
        downsampleStatsTo: 'techrepDownsampleStatsTo',
        downsamplePairs: 'techrepDownsamplePairs',
        doMerge: 'mergeTechrepToBiorep',
        doDedup: 'techrepDedup',
        mergeGroupIdentifiers: ['condition', 'biorep', 'aggregateProfileName'],
        levelFilter: {it.condition && it.biorep && it.techrep}
    ]

    samples | branch {
        techrep: isTechrep(it)
        other: true
    } | set{sampleType}

    CreateAggregatePairProfiles(sampleType.techrep, aggregateParams)
    | DownsamplePairs
    | MergePairs
    | DeduplicatePairs
    | concat(sampleType.other)
    | set{result}


    emit:
    result
}

workflow AggregateBioreps {
    take:
    samples

    main:

    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'max',
        cisStrata: 'biorepCisStrata',
        readConjuncts: 'biorepReadConjuncts',
        downsampleToMeanDistribution: 'biorepDownsampleToMeanDistribution',
        downsampleToSize: 'biorepDownsampleToSize',
        downsampleStatsFrom: 'biorepDownsampleStatsFrom',
        downsampleStatsTo: 'biorepDownsampleStatsTo',
        downsamplePairs: 'biorepDownsamplePairs',
        doMerge: 'mergeBiorepToCondition',
        doDedup: 'biorepDedup',
        mergeGroupIdentifiers: ['condition', 'aggregateProfileName'],
        levelFilter: {it.condition && it.biorep && !it.techrep}
    ]

    samples | branch {
        biorep: isBiorep(it)
        other: true
    } | set{sampleType}

    CreateAggregatePairProfiles(sampleType.biorep, aggregateParams)
    | DownsamplePairs
    | MergePairs
    | DeduplicatePairs
    | concat(sampleType.other)
    | set{result}

    emit:
    result
}

workflow AggregateConditions {
    take:
    samples

    main:

    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'max',
        cisStrata: 'conditionCisStrata',
        readConjuncts: 'conditionReadConjuncts',
        downsampleToMeanDistribution: 'conditionDownsampleToMeanDistribution',
        downsampleToSize: 'conditionDownsampleToSize',
        downsampleStatsFrom: 'conditionDownsampleStatsFrom',
        downsampleStatsTo: 'conditionDownsampleStatsTo',
        downsamplePairs: 'conditionDownsamplePairs',
        doMerge: 'mergeCondition',
        doDedup: 'conditionDedup',
        mergeGroupIdentifiers: ['aggregateProfileName'],
        levelFilter: {it.condition && !it.biorep && !it.techrep}
    ]

    samples | branch {
        condition: isCondition(it)
        other: true
    } | set{sampleType}

    CreateAggregatePairProfiles(sampleType.condition, aggregateParams)
    | DownsamplePairs
    | MergePairs
    | DeduplicatePairs
    | concat(sampleType.other)
    | set{result}

    emit:
    result
}
