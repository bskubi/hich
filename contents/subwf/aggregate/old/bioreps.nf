include {CreateAggregationPlans} from './aggregation_plans.nf'

workflow AggregateBioreps {
    take:
    samples

    main:

    // Separate the biorep-level samples where we are to aggregate the bioreps and that have pairs
    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'biorepDedupMethod',
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

    // Separate the bioreps-level samples where we are to aggregate the bioreps and that have pairs
    samples | branch {
        biorep: !skip("aggregate") && !skip("aggregateBioreps") && isBiorep(it) && (it.pairs || it.latestPairs)
        other: true
    } | set{sampleType}

    // Downsample, merge and deduplicate as necessary
    CreateAggregationPlans(sampleType.bioreps, aggregateParams) | view

    // Early stopping
    samples = emptyOnLastStep("aggregateBioreps", samples)

    emit:
    samples
}