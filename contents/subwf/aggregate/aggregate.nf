include {emptyOnLastStep} from '../util/workflowControl.nf'
include {processMerge} from './processMerge.nf'
include {Merge as MergeTechrepsToBioreps; Merge as MergeBiorepsToConditions} from './merge.nf'
include {Deduplicate} from './deduplicate.nf'
include {Split} from './split.nf'
include {columnsToRows} from '../util/rowsCols.nf'
include {pack} from '../util/join.nf'
include {makeID} from '../util/samples.nf'


workflow Aggregate {
    take:
    samples

    main:

    // Merge samples
    samples
    | branch {
        techreps: it.aggregateLevel == "techrep"
        bioreps: it.aggregateLevel == "biorep"
        conditions: it.aggregateLevel == "condition"
    }
    | set {sampleLevels}

    processMerge(
        sampleLevels.techreps, 
        ["cell", "biorep", "condition", "aggregationPlanName"], 
        MergeTechrepsToBioreps,
        "biorep",
        sampleLevels.bioreps
        )
    | set{bioreps}

    samples = emptyOnLastStep("mergeTechrepsToBioreps", samples)

    processMerge(
        bioreps, 
        ["cell", "condition", "aggregationPlanName"], 
        MergeBiorepsToConditions,
        "condition",
        sampleLevels.conditions
        )
    | set{conditions}

    samples = emptyOnLastStep("mergeBiorepsToConditions", samples)
    samples = emptyOnLastStep("merge", samples)

    sampleLevels.techreps
    | concat(bioreps)
    | concat(conditions)
    | set{samples}

    // Deduplicate samples
    samples
    | map{tuple(it.id, it.latestPairs, it.dedupSingleCell, it.dedupMaxMismatch, it.dedupMethod, it.pairtoolsDedupParams)}
    | Deduplicate
    | map{[id:it[0], dedupPairs:it[1], latest:it[1], latestPairs:it[1]]}
    | set {deduplicated}

    pack(samples, deduplicated)
    | set {samples}

    samples = emptyOnLastStep("deduplicate", samples)

    // Split samples based on cell column
    samples
    | map{tuple(it.id, it.latestPairs)}
    | Split
    | map{[id: [it[0]]*it[1].size(), splitPairs: it[1], latestPairs: it[1]]}
    | columnsToRows
    | set{splitSamples}

    pack(splitSamples, samples)
    | map{
        cell = it.splitPairs =~ /.*\.cell=([^.]+)\..*/

        it += [cell: cell.matches() ? cell[0][1] : null ]
        it += [id: makeID(it)]
    }
    | set{samplesFromSplit}

    samples
    | concat(samplesFromSplit)
    | set{samples}

    samples = emptyOnLastStep("split", samples)

    samples = emptyOnLastStep("aggregate", samples)

    emit:
    samples
}