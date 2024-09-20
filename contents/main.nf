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

/*
Splitting the aggregate file into individual data files is solved cleanly via hich organize

Building and loading in the data files to Hich is not.

There are several conceptual questions:
    1. When in the workflow to do the splitting and loading in
        I think as soon as possible after alignment
            This takes advantage of paralellism to the maximum extent, avoids unnecessary deduplication,
            avoids unnecessary carryover of sam/bam information
    2. Whether the techrep, biorep and condition fields are sufficient or whether we need a "cell" or "cell type" field as well
       for single-cell data
        I think having an additional 1-2 fields for single-cell is appropriate (?)
            Imagine we have two drug treatments with a single-cell assay (conditions)
            We do multiple specimen samples for each treatment (bioreps)
            We sequence in multiple batches (techreps)
            Then when we split into single cells, we would have to discard one of these pieces of information
            The other alternative would be to think of each cell or cell type as a biorep
            This involves erasing the difference between the specimen samples, which seems inappropriate

            Can we use the cell during pseudobulking? I think so -- just treat it as a "bigger cell"
            A cell, then, is just an extra tag that we can potentially merge on and forms part of the id IF given.
            Doesn't have to be given for bulk data.
            It can be directly produced from a split or a merge.
            A pseudobulk cell 


    3. How to obtain parameters from the combination of the aggregate file, generated sample file, and nextflow params


*/

workflow {
    LoadSampleFile              // Setup workflow inputs
        | Setup
        | FastqHead

        | Align                 // Align .fastq -> .bam
        | Parse                 // .bam -> .pairs and read-level filters
        | IngestPairs
        | TagFragments
        | Select

        //| DownsampleTechreps    // Condition/biorep/techrep coverage control and merge
        | TechrepsToBioreps
        | Deduplicate
        //| DownsampleBioreps
        | BiorepsToConditions
        //| DownsampleConditions
        
        | HicMatrix             // Create contact matrices
        | McoolMatrix

        | Hicrep                // Call features
        | CompartmentScore
        | DifferentialLoops
        | InsulationScore
        | set{samples}

    samples = emptyOnLastStep("End", samples)
}

