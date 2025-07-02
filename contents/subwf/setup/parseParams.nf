include {emptyOnLastStep} from '../util/cli.nf'
include {withLog; stubLog} from '../util/logs.nf'
include {ParseFastqInterleaved} from './parseParams/parseFastqInterleaved.nf'
include {ParseFastqPairs} from './parseParams/parseFastqPairs.nf'
include {ParseSampleFile} from './parseParams/parseSampleFile.nf'
include {ParseSamples} from './parseParams/parseSamples.nf'
include {ParseSamplesFromSRA} from './parseParams/parseSamplesFromSRA.nf'
 
workflow ParseParams {
    main:
    channel.empty()
        | ParseSampleFile
        | ParseFastqInterleaved
        | ParseFastqPairs
        | ParseSamples
        | ParseSamplesFromSRA
        | set {samples}

    samples = emptyOnLastStep("parseParams", samples)

    emit:
    samples
}