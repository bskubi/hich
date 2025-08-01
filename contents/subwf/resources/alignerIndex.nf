include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {PrepIndex} from './helper/prepIndex.nf'
include {BwaMemIndex} from './processes/bwaMemIndex.nf'
include {BwaMem2Index} from './processes/bwaMem2Index.nf'
include {BwamethMem2Index} from './processes/bwamethMem2Index.nf'
include {BwamethIndex} from './processes/bwamethIndex.nf'

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
                    [genomeReference: file(genomeReference), alignerIndexDir: alignerIndexDir, alignerIndexPrefix: alignerIndexPrefix, aligner: "bwa-mem2"]}
            | set{resultBwaMem2Index}

        samples
            | filter{it.aligner == "bwa"}
            | PrepIndex
            | BwaMemIndex
            | map{genomeReference, alignerIndexDir, alignerIndexPrefix, prefix_ann, prefix_amb, prefix_pac, prefix_bwt, prefix_sa ->
                    [genomeReference: file(genomeReference), alignerIndexDir: alignerIndexDir, alignerIndexPrefix: alignerIndexPrefix, aligner: "bwa"]}
            | set{resultBwaMemIndex}

        samples
            | filter{it.aligner == "bwameth-mem2"}
            | PrepIndex
            | BwamethMem2Index
            | map{genomeReference, alignerIndexDir, alignerIndexPrefix,
                c2t, c2t_amb, c2t_ann, c2t_bwt_2bit_64, c2t_pac, c2t_0123 ->
                    [genomeReference: file(genomeReference), alignerIndexDir: alignerIndexDir, alignerIndexPrefix: alignerIndexPrefix, aligner: "bwameth-mem2"]}
            | set{resultBwamethMem2Index}

        samples
            | filter{it.aligner == "bwameth"}
            | PrepIndex
            | BwamethIndex
            | map{genomeReference, alignerIndexDir, alignerIndexPrefix,
                c2t, c2t_amb, c2t_ann, c2t_bwt, c2t_pac, c2t_sa ->
                    [genomeReference: file(genomeReference), alignerIndexDir: alignerIndexDir, alignerIndexPrefix: alignerIndexPrefix, aligner: "bwameth"]}
            | set{resultBwamethIndex}

            joinOn = ["genomeReference", "aligner"]

        // TODO: Skip samples with an aligner index and prefix already given

        keyUpdate(samples, resultBwaMem2Index, joinOn) | set{samples}
        keyUpdate(samples, resultBwaMemIndex, joinOn) | set{samples}
        keyUpdate(samples, resultBwamethMem2Index, joinOn) | set{samples}
        keyUpdate(samples, resultBwamethIndex, joinOn) | set{samples}
    }


    samples = emptyOnLastStep("alignerIndex", samples)


    emit:
        samples

}

