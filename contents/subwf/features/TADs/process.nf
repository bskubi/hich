include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process TADS {
    publishDir "results/tads",
               mode: params.general.publish.mode
    label 'features'
    container params.general.hicexplorerContainer
    tag "$id"

    input:
    tuple val(id), path(mcool), val(tads_opts)

    output:
    tuple val(id), path(boundariesBed), path(boundariesGff), path(domainsBed), path(scoreBedgraph), path(tadScoreBm)

    shell:
    (cmd, logMap, boundariesBed, boundariesGff, domainsBed, scoreBedgraph, tadScoreBm) = buildCmd(id, mcool, tads_opts, task.cpus)

    
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, boundariesBed, boundariesGff, domainsBed, scoreBedgraph, tadScoreBm) = buildCmd(id, mcool, tads_opts, task.cpus)
    stub = "touch '${output.boundariesBed}' '${output.boundariesGff}' '${output.domainsBed}' '${output.scoreBedgraph}' '${output.tadScoreBm}'"
    stubLog(stub, cmd, logMap)
}