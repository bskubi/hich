include {ParseParams} from './subwf/setup/parseParams.nf'
include {Setup} from './subwf/setup/setup.nf'

include {Align} from './subwf/reads/Align/workflow.nf'
include {ParseToPairs} from './subwf/reads/ParseToPairs/workflow.nf'
include {TagRestrictionFragments} from './subwf/reads/TagRestrictionFragments/workflow.nf'
include {IngestPairs} from './subwf/reads/IngestPairs/workflow.nf'
include {SelectPairs} from './subwf/reads/SelectPairs/workflow.nf'

include {AggregatePairs} from './subwf/aggregate/workflow.nf'

include {CreateMatrix} from './subwf/matrix/workflow.nf'
include {Hicrep} from './subwf/features/hicrep.nf'
include {CompartmentScores} from './subwf/features/CompartmentScores/workflow.nf'
include {InsulationScores} from './subwf/features/insulationScores.nf'
include {TADs} from './subwf/features/tads.nf'
include {Loops} from './subwf/features/Loops/workflow.nf'
include {DifferentialLoops} from './subwf/features/features/DifferentialLoops/workflow.nf'

include {emptyOnLastStep; skip} from './subwf/util/cli.nf'

workflow HichWorkflow {
    
    ParseParams
        | Setup

        | Align
        | SambamToPairs
        | IngestPairs
        | TagRestrictionFragments
        | SelectPairs

        | AggregatePairs
        | CreateMatrix

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

