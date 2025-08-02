include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process TAG_RESTRICTION_FRAGMENTS {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null},
               mode: params.general.publish.mode

    label 'pairs'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer

    input:
    tuple val(id), path(pairs), path(fragmentIndex)

    output:
    tuple val(id), path(output)

    shell:
    (cmd, logMap, output) = buildCmd(id, pairs, fragmentIndex)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmd(id, pairs, fragmentIndex)
    stub = "touch '${id}_fragtag.pairs.gz'"
    stubLog(stub, cmd, logMap)
}