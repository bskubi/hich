include {createCompositeStrategy; filterSamplesByStrategy} from '../../util/analysisPlans.nf'
include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {MCOOL_MATRIX} from './process.nf'

workflow McoolMatrix {
    take:
    samples
    
    main:
    myName = "mcoolMatrix"
    if (!skip(myName)) {

        samples
            | filter{it.makeMcoolFileFormat && it.latestPairs && !it.mcool}
            | map{tuple(it.id, it.latestPairs, it.chromsizes, it.assembly, it.matrix_opts)}
            | MCOOL_MATRIX
            | map{
                id, mcool, matrix_opts -> 
                [id: id, matrix_opts: matrix_opts, mcool: mcool, latestMatrix: mcool]
            }
            | set{result}
        
        keyUpdate(samples, result, ["id", "matrix_opts"])
            | set{samples}
    }
    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}