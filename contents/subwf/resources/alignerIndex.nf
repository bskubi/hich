include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {isExistingFile} from '../util/files.nf'
include {withLog; stubLog} from '../util/logs.nf'
include {BSBoltIndex} from './processes/bsboltIndex.nf'
include {BwaMemIndex} from './processes/bwaMemIndex.nf'
include {BwaMem2Index} from './processes/bwaMem2Index.nf'

workflow PrepIndex {
    take:
        samples
    main:
        samples
            | filter {it.datatype == "fastq"}
            | filter {!isExistingFile(it.alignerIndexDir)}
            | map{it.alignerIndexPrefix = it.alignerIndexPrefix ?: it.assembly; it}
            | map{tuple(it.genomeReference, it.genomeReference, it.alignerIndexPrefix)}
            | unique
            | set{result}
    emit:
        result
}

workflow AlignerIndex {

    take:
        samples
    
    main:

    if (!skip("alignerIndex")) {

    samples
        | filter{it.aligner == "bwa-mem2"}
        | PrepIndex
        | BwaMem2Index
        | map{genomeReference, alignerIndexDir, alignerIndexPrefix, prefix_0123, prefix_amb, prefix_ann, prefix_bwt_2bit_64, prefix_pac ->
                [genomeReference: file(genomeReference), alignerIndexDir: alignerIndexDir, alignerIndexPrefix: alignerIndexPrefix]}
        | set{resultBwaMem2Index}

    samples
        | filter{it.aligner == "bwa"}
        | PrepIndex
        | BwaMemIndex
        | map{genomeReference, alignerIndexDir, alignerIndexPrefix, prefix_ann, prefix_amb, prefix_pac, prefix_bwt, prefix_sa ->
                [genomeReference: file(genomeReference), alignerIndexDir: alignerIndexDir, alignerIndexPrefix: alignerIndexPrefix]}
        | set{resultBwaMemIndex}

    // TODO add bwameth indexing

    samples
        | filter {it.datatype == "fastq" && it.aligner == "bsbolt"}
        | PrepIndex
        | BSBoltIndex
        | map{genomeReference, alignerIndexDir, alignerIndexPrefix,
              prefix_fa, prefix_ann, prefix_amb, prefix_opac, prefix_pac, prefix_bwt, prefix_sa ->
                [genomeReference: file(genomeReference), alignerIndexDir: alignerIndexDir, alignerIndexPrefix: alignerIndexPrefix]}
        | set{resultBSBoltIndex}

        keyUpdate(samples, resultBwaMem2Index, "genomeReference") | set{samples}
        keyUpdate(samples, resultBwaMemIndex, "genomeReference") | set{samples}
        keyUpdate(samples, resultBSBoltIndex, "genomeReference") | set{samples}
    }


    samples = emptyOnLastStep("alignerIndex", samples)

    emit:
        samples

}

