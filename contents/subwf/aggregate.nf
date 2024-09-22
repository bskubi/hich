include {updateChannel; coalesce; rows; columns; groupHashMap; isTechrep; constructIdentifier} from './extraops.nf'

process HichDownsampleStatsFrom {
    input:
    tuple val(id), val(conjuncts), val(cisStrata)

    output:
    tuple val(id), path(stats)

    stub:
    stats = "${id}.stats.tsv"
    "touch ${stats}"
}

process HichDownsampleStatsTo {
    input:
    tuple val(ids), path(stats), val(toSizes), val(outliers)

    output:
    tuple val(ids), path(targetStats)

    shell:
    targetStats = ids.collect{"${it}.target.stats.tsv"}
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
*/

workflow AggregateTechreps {
take:
samples

main:

samples | filter{isTechrep(it) && it.aggregateProfile == null} | set{rawTechreps}


channel.empty() | set{downsampleTechreps}
params.comparisonSets.aggregateProfiles.each {
    profile, profileParams ->

    
    rawTechreps
    | map{it + [aggregateProfile:profile] + profileParams}
    | map{it + [techrepToMeanDistribution: it.subMap(it.techrepToMeanDistribution) ?: null] + [id: "${it.id}_${it.aggregateProfile}"]}
    | concat(downsampleTechreps)
    | set{downsampleTechreps}
}

downsampleTechreps
    | map{tuple(it.id, it.techrepReadConjuncts, it.techrepCisStrata)}
    | HichDownsampleStatsFrom
    | map{[id:it[0], techrepDownsampleStatsFrom:it[1]]}
    | set{downsampleStatsFromResult}
updateChannel(downsampleTechreps, downsampleStatsFromResult) | set {fromStatsCalculated}

groupHashMap(fromStatsCalculated, 'techrepToMeanDistribution')
    | map{columns(it)}
    | map{tuple(it.id, it.techrepDownsampleStatsFrom, it.techrepToSize, it.outlier)}
    | HichDownsampleStatsTo
    | map{[id:it[0], techrepDownsampleStatsTo:it[1]]}
    | map{rows(it)}
    | flatten
    | set {downsampleStatsToResult}
updateChannel(fromStatsCalculated, downsampleStatsToResult) | set{toStatsCalculated}

toStatsCalculated
    | map{tuple(it.id, it.techrepReadConjuncts, it.techrepCisStrata, it.techrepDownsampleStatsFrom, it.techrepDownsampleStatsTo)}
    | HichDownsamplePairs
    | map{[id:it[0], techrepDownsampledPairs:it[1], latest:it[1], latestPairs:it[1]]}
    | set{downsampledResult}
updateChannel(toStatsCalculated, downsampledResult) | set {downsampled}

downsampled | branch {
    yes: it.includeInMerge != false
    no: true}
| set{merge}



groupHashMap(merge.yes, ['condition', 'biorep', 'aggregateProfile'])
    | map{coalesce(columns(it))}

    | map{tuple(constructIdentifier(coalesce(it, "_drop")), it.latestPairs)}
    | PairtoolsMerge
    | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
    | set{biorepsFromMerge}


groupHashMap(merge.yes, ['condition', 'biorep', 'aggregateProfile'])
    | map{coalesce(columns(it), '_drop')}
    | map{it += [id:constructIdentifier(it)]}
    | set{inheritedBiorepsAttributes}
updateChannel(inheritedBiorepsAttributes, biorepsFromMerge) | set{biorepsFromMerge}

rawTechreps
    | concat(merge.yes, merge.no)
    | branch {
    yes: it.deduplicate
    no: true
} | set{deduplicate}

deduplicate.yes
    | map{tuple(it.id, it.latestPairs, it.singleCell, it.dedupMaxMismatch, it.dedupMethod)}
    | PairtoolsDedup
    | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
    | set{dedupResult}
updateChannel(deduplicate.yes, dedupResult) | concat(deduplicate.no, biorepsFromMerge) | set{samples}

emit:
samples
}