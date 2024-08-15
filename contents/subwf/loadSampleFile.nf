include {emptyOnLastStep} from './extraops.nf'

workflow LoadSampleFile {
    main:
    samples = channel.fromPath(params.sampleFile, checkIfExists: true)
        | splitCsv(header: true, sep: params.sampleFileSep)
        
    samples = emptyOnLastStep("loadSampleFile") ?: samples

    emit:
    samples
}