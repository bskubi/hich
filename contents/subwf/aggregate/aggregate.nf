
include {Merge as MergeTechrepToBiorep; Merge as MergeBiorepToCondition} from './merge.nf'
include {LabelAggregationPlans} from './labelAggregationPlans.nf'
include {Deduplicate} from './deduplicate.nf'
include {Split} from './split.nf'
include {QCPairs as DedupQCPairs; QCPairs as AggregateQCPairs} from '../reads/qcPairs.nf'
include {emptyOnLastStep; skip} from '../util/cli.nf'
include {columnsToRows} from '../util/reshape.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {makeID} from '../util/samples.nf'


workflow Aggregate {
    take:
    samples

    main:

    if (!skip("aggregate")) {
        samples
            | LabelAggregationPlans
            | set{samples}

        if (!skip("mergeTechreps") && !skip("merge")) {
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
        if (!skip("dedup")) {
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
            | Deduplicate
            | map{[id:it[0], dedupPairs:it[1], latest:it[1], latestPairs:it[1]]}
            | set {deduplicated}

            keyUpdate(dedup.yes, deduplicated, "id")
            | concat(dedup.no)
            | set {samples}

            if ("dedup" in params.general.get("qcAfter")) {
                DedupQCPairs(samples, ["dedupPairs"], "dedup")
            }
        }

        samples = emptyOnLastStep("dedup", samples)

        // Merge bioreps to conditions

        if (!skip("mergeBiorepToCondition") && !skip("merge")) {
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

        samples = emptyOnLastStep("mergeBiorepToCondition", samples)

        // Split into multiple samples (i.e. on cell ID)

        if (!skip("split")) {
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
                | Split
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

        if ("aggregate" in params.general.get("qcAfter")) {
            AggregateQCPairs(samples, ["latestPairs"], "aggregate")
        }

    }

    samples = emptyOnLastStep("aggregate", samples)

    emit:
    samples
}