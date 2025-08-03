
include {MergePairs as MergeTechrepToBiorep; MergePairs as MergeBiorepToCondition} from './MergePairs/workflow.nf'
include {LabelAggregationPlans} from './LabelAggregationPlans/workflow.nf'
include {DEDUP_PAIRS} from './DedupPairs/process.nf'
include {SPLIT_PAIRS} from './SplitPairs/process.nf'
include {QCPairs as DedupQCPairs; QCPairs as AggregateQCPairs} from '../reads/QCPairs/workflow.nf'
include {emptyOnLastStep; skip} from '../util/cli.nf'
include {columnsToRows} from '../util/reshape.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {makeID} from '../util/samples.nf'


workflow Aggregate {
    take:
    samples

    main:

    if (!skip("Aggregate")) {
        samples
            | LabelAggregationPlans
            | set{samples}

        if (!skip("MergeTechreps") && !skip("merge")) {
            samples
                | branch {
                    yes: it.aggregateLevel == "techrep" && !it.skipMerge && !it.skipTechrepMerge && it.mergeTechrepToBiorep
                    no: true
                }
                | set {mergeTechreps}


            // Merge techreps to bioreps
            MergeTechrepToBiorep(
                mergeTechreps.yes, 
                ["cell", "biorep", "condition", "aggregationPlanName"],
                "biorep"
                )
            | LabelAggregationPlans
            | concat(samples)
            | set {samples}

        }

        samples = emptyOnLastStep("mergeTechrepToBiorep", samples)

        // Deduplicate techreps, bioreps, and input conditions
        if (!skip("DedupPairs")) {
            samples
            | branch {
                yes: (
                    it.dedupMethod 
                    || it.dedupMaxMismatch 
                    || it.dedupSingleCell 
                    || it.dedup
                    || (it.aggregateLevel == "techrep" && it.dedupTechreps)
                    || (it.aggregateLevel == "biorep" && it.dedupBioreps)
                    || (it.aggregateLevel == "condition" && it.dedupConditions)
                )
                no: true
            }
            | set{dedup}

            dedup.yes
            | map{tuple(it.id, it.latestPairs, it.dedupSingleCell, it.dedupMaxMismatch, it.dedupMethod, it.pairtoolsDedupParams)}
            | DEDUP_PAIRS
            | map{[id:it[0], dedupPairs:it[1], latest:it[1], latestPairs:it[1]]}
            | set {deduplicated}

            keyUpdate(dedup.yes, deduplicated, "id")
            | concat(dedup.no)
            | set {samples}

            if ("DedupPairs" in params.general.get("qcAfter")) {
                DedupQCPairs(samples, ["dedupPairs"], "DedupPairs")
            }
        }

        samples = emptyOnLastStep("dedup", samples)

        // Merge bioreps to conditions

        if (!skip("MergeBiorepToCondition") && !skip("merge")) {
            samples
                | branch {
                    yes: it.aggregateLevel == "biorep" && !it.skipMerge && !it.skipBiorepMerge && it.mergeBiorepToCondition
                    no: true
                }
                | set {mergeBioreps}

            MergeBiorepToCondition(
                mergeBioreps.yes, 
                ["cell", "condition", "aggregationPlanName"],
                "condition"
                )
                | LabelAggregationPlans
                | concat(samples)
                | set{samples}

        }

        samples = emptyOnLastStep("MergeBiorepToCondition", samples)

        // Split into multiple samples (i.e. on cell ID)

        if (!skip("SplitPairs")) {
            samples
                | branch {
                    yes: (
                        it.split 
                        || it.splitColumns
                        || it.splitSQL
                        || (it.aggregateLevel == "techrep" && it.splitTechreps)
                        || (it.aggregateLevel == "biorep" && it.splitBioreps)
                        || (it.aggregateLevel == "condition" && it.splitConditions)
                    )
                    no: true
                }
                | set {split}

            split.yes
                | map{tuple(it.id, it.latestPairs, it.splitColumns, it.splitSQL)}
                | SPLIT_PAIRS
                | map{[id: [it[0]]*it[1].size(), splitPairs: it[1], latestPairs: it[1]]}
                | columnsToRows
                | set{splitSamples}

            keyUpdate(splitSamples, split.yes, "id")
                | map{
                    // Extract cell from filename
                    cell = it.splitPairs =~ /.*\.cell=([^.]+)\..*/

                    it += [cell: cell.matches() ? cell[0][1] : null ]
                    it += [id: makeID(it, false)]
                }
                | LabelAggregationPlans
                | concat(samples)
                | set{samples}
        }

        if ("Aggregate" in params.general.get("qcAfter")) {
            AggregateQCPairs(samples, ["latestPairs"], "Aggregate")
        }

    }

    samples = emptyOnLastStep("Aggregate", samples)

    emit:
    samples
}