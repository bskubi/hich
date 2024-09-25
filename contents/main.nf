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
include {AggregateTechreps; AggregateBioreps; AggregateConditions} from './subwf/aggregate.nf'


workflow {
    
    LoadSampleFile              // Setup workflow inputs
        | Setup
        | FastqHead

        | Align                 // Align .fastq -> .bam
        | Parse                 // .bam -> .pairs and read-level filters
        | IngestPairs
        | TagFragments
        | Select

        | AggregateTechreps
        | AggregateBioreps
        | AggregateConditions
        
        | HicMatrix             // Create contact matrices
        | McoolMatrix

        | Hicrep                // Call features
        | CompartmentScore
        | DifferentialLoops
        | InsulationScore

        | set{samples}

    samples = emptyOnLastStep("End", samples)
}

