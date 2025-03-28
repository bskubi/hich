include {CreateAggregationPlans} from './aggregation_plans.nf'

workflow AggregateConditions {
    take:
    samples

    main:

    // Separate the condition-level samples where we are to aggregate the bioreps and that have pairs
    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'conditionDedupMethod',
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

    // Separate the condition-level samples where we are to aggregate the bioreps and that have pairs
    samples | branch {
        condition: !skip("aggregate") && !skip("aggregateConditions") && isCondition(it) && (it.pairs || it.latestPairs)
        other: true
    } | set{sampleType}
    
    // Downsample, merge and deduplicate as necessary
    CreateAggregationPlans(sampleType.condition, aggregateParams) | view

    // Early stopping
    samples = emptyOnLastStep("aggregateConditions", samples)

    emit:
    samples
}
