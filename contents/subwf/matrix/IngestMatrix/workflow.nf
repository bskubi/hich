include {emptyOnLastStep; skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {HIC_TO_MCOOL; MCOOL_TO_HIC} from './process.nf'


workflow IngestMatrix {
    take:
        samples

    main:

    if (!skip("IngestMatrix")) {
        samples
            | filter{it.datatype == "matrix" && it.makeHicFileFormat && it.mcool && !it.hic}
            | map{tuple(it.id, it.mcool)}
            | MCOOL_TO_HIC
            | map{
                id, hicFile ->
                [id:id, hic:hicFile, latest:hicFile, latestMatrix:hicFile]}
            | set{mcoolToHic}
        
        samples
            | filter{it.datatype == "matrix" && it.makeMcoolFileFormat && it.hic && !it.mcool}
            | map{tuple(it.id, it.hic)}
            | HIC_TO_MCOOL
            | map{
                id, mcoolFile ->
                [id:id, mcool:mcoolFile, latest:mcoolFile, latestMatrix:mcoolFile]}
            | set{hicToMcool}

        keyUpdate(samples, mcoolToHic, "id") | set{samples}
        keyUpdate(samples, hicToMcool, "id") | set{samples}
    }

    


    samples = emptyOnLastStep("IngestMatrix", samples)

    emit:
        samples
}