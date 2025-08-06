include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {JUICER_TOOLS_PRE} from './process.nf'

workflow HicMatrix {
    take:
    samples
    
    main:
    myName = "HicMatrix"
    if (!skip(myName)) {
        samples
            | filter{it.makeHicFileFormat && it.latestPairs && !it.hic}
            | map{tuple(it.id, it.latestPairs, it.chromsizes, it.matrix_opts, it.minMapq)}
            | JUICER_TOOLS_PRE
            | map{
                id, matrix_opts, hic -> 
                [id: id, matrix_opts: matrix_opts, hic: hic, latestMatrix: hic]
            }
            | set{result}
            
            keyUpdate(samples, result, ["id", "matrix_opts"])
                | set{samples}
    }


    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}
