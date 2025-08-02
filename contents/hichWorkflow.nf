include {ParseParams} from './subwf/setup/parseParams.nf'
include {Setup} from './subwf/setup/setup.nf'

include {Align} from './subwf/reads/align.nf'
include {Parse} from './subwf/reads/parse.nf'
include {TagRestrictionFragments} from './subwf/reads/TagRestrictionFragments/workflow.nf'
include {IngestPairs} from './subwf/reads/ingestPairs.nf'
include {Select} from './subwf/reads/Select/workflow.nf'

include {LabelAggregationPlans} from './subwf/aggregate/labelAggregationPlans.nf'
include {Aggregate} from './subwf/aggregate/aggregate.nf'

include {LabelMatrixPlans} from './subwf/matrix/labelMatrixPlans.nf'
include {HicMatrix} from './subwf/matrix/hicMatrix.nf'
include {McoolMatrix} from './subwf/matrix/mcoolMatrix.nf'
include {IngestMatrix} from './subwf/matrix/ingestMatrix.nf'
include {Hicrep} from './subwf/features/hicrep.nf'
include {CompartmentScores} from './subwf/features/compartmentScores.nf'
include {InsulationScores} from './subwf/features/insulationScores.nf'
include {TADs} from './subwf/features/tads.nf'
include {Loops} from './subwf/features/loops.nf'
include {DifferentialLoops} from './subwf/features/differentialLoops.nf'

include {emptyOnLastStep; skip} from './subwf/util/cli.nf'

workflow HichWorkflow {
    
    ParseParams
        | Setup

        | Align
        | Parse
        | IngestPairs
        | TagRestrictionFragments
        | Select

        | Aggregate

        | LabelMatrixPlans
        | IngestMatrix
        | HicMatrix
        | McoolMatrix

        | Hicrep
        | CompartmentScores
        | InsulationScores
        | TADs
        | Loops
        | DifferentialLoops

        | set{samples}

    emit:
    samples
}

