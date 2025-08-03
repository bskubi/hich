include {LabelAggregationPlans} from '../LabelAggregationPlans/workflow.nf'
include {SPLIT_PAIRS} from './process.nf'
include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {columnsToRows} from '../../util/reshape.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {makeID} from '../../util/samples.nf'

workflow SplitPairs {
    take:
    samples

    main:
    myName = "SplitPairs"

    if (!skip(myName)) {
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

    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}