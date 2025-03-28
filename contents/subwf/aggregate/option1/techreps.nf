include {CreateAggregationPlans} from './aggregation_plans.nf'
include {emptyOnLastStep; skip} from '../util/workflow_control.nf'
include {isTechrep} from '../util/samples.nf'

workflow AggregateTechreps {
    take:
    samples

    main:

    // Create techreps-specific versions of the generic parameters
    aggregationParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'techrepDedupMethod',
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

    // Separate the techrep-level samples where we are to aggregate the techreps and that have pairs
    samples
    | branch {
        techrep: !skip("aggregate") && !skip("aggregateTechreps") && isTechrep(it) && (it.pairs || it.latestPairs)
        other: true
    } 
    | set{sampleType}

    // Downsample, merge and deduplicate as necessary
    CreateAggregationPlans(sampleType.techrep, aggregationParams)
    | concat(sampleType.other)
    | set{samples}

    // Early stopping
    samples = emptyOnLastStep("aggregateTechreps", samples)

    emit:
    samples
}