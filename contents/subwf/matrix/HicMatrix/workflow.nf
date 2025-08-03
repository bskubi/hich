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
            | map{tuple(it.id, it.matrixPlanName, it.latestPairs, it.chromsizes, it.matrix, it.juicerToolsPreParams, it.minMapq)}
            | JUICER_TOOLS_PRE
            | map{
                id, matrixPlanName, hic -> 
                [id: id, matrixPlanName: matrixPlanName, hic: hic, latestMatrix: hic]
            }
            | set{result}
            
            keyUpdate(samples, result, ["id", "matrixPlanName"])
                | set{samples}
    }


    samples = emptyOnLastStep(myName, samples)

    emit:
    samples
}
