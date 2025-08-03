include {skip} from '../../util/cli.nf'
include {keyUpdate} from '../../util/keyUpdate.nf'
include {PAIRTOOLS_STATS; MULTIQC_PAIRTOOLS} from './process.nf'

workflow QCPairs {
    take:
    samples
    pairs
    reportName
    
    main:

    if (!skip("qcpairs")) {
        samples
            | map{
                sample ->
                
                pairs.collect{
                    pairfile ->
                    hasPairFile = sample.containsKey(pairfile)
                    pairFileIsPath = sample.get(pairfile).getClass() == sun.nio.fs.UnixPath
                    statsKey = pairfile + "Stats"
                    noStatsFile = !sample.containsKey(statsKey)
                    computeStats = hasPairFile && pairFileIsPath && noStatsFile
                    statsInput = [sample.id, pairfile, statsKey, sample[pairfile]]
                    computeStats ? statsInput : null
                }.findAll{it != null}
            }
            | collect
            | flatMap
            | PAIRTOOLS_STATS
            | map{it[3]}
            | collect
            | map{paths -> [reportName, paths]}
            | MULTIQC_PAIRTOOLS
    }



    emit:
    samples

}