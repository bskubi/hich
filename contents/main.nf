include {LoadSampleFile} from './subwf/loadSampleFile.nf'
include {Setup} from './subwf/setup.nf'
include {Align} from './subwf/align.nf'
include {FastqHead} from './subwf/fastqHead.nf'
include {Parse} from './subwf/parse.nf'
include {TagFragments} from './subwf/tagFragments.nf'
include {TechrepsToBioreps; BiorepsToConditions} from './subwf/mergePairs.nf'
include {IngestPairs} from './subwf/ingestPairs.nf'
include {Deduplicate} from './subwf/deduplicate.nf'
include {Select} from './subwf/select.nf'
include {HicMatrix} from './subwf/hicMatrix.nf'
include {McoolMatrix} from './subwf/mcoolMatrix.nf'
include {Hicrep} from './subwf/hicrep.nf'
include {CompartmentScore} from './subwf/compartmentScore.nf'
include {DifferentialLoops} from './subwf/differentialLoops.nf'
include {InsulationScore} from './subwf/insulationScore.nf'
include {emptyOnLastStep} from './subwf/extraops.nf'

workflow {
    LoadSampleFile
        | Setup
        | FastqHead
        | Align
        | Parse
        | IngestPairs
        | TagFragments
        | TechrepsToBioreps
        | Deduplicate
        | BiorepsToConditions
        | Select
        | HicMatrix
        | McoolMatrix
        | Hicrep
        | CompartmentScore
        | DifferentialLoops
        | InsulationScore
        | set{samples}
    samples = emptyOnLastStep("End", samples)
}

