include {withLog; stubLog} from '../../util/logs.nf'
include {buildCmd} from './functions.nf'

process COMPARTMENT_SCORES {
    publishDir "results/compartments",
               mode: params.general.publish.mode

    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer
    label 'features'
    tag "$id"


    input:
    tuple val(id), path(genomeReference), path(matrix), val(resolution), val(eigsCisParams), val(nEigs)

    output:
    tuple val(id), path("*.cis.bw"), path("*.cis.vecs.tsv"), path("*.cis.lam.txt"), path("*.phase.bed")

    shell:
    (cmd, logMap, output) = buildCmd(id, genomeReference, matrix, resolution, eigsCisParams, nEigs)
    withLog(cmd, logMap)

    stub:
    (cmd, logMap, output) = buildCmd(id, genomeReference, matrix, resolution, eigsCisParams, nEigs)
    stub = "touch '${id}_compartments.bw' '${id}.cis.bw' '${id}.cis.vecs.tsv' '${id}.cis.lam.txt' '${id}.phase.bed'"
    stubLog(stub, cmd, logMap)
}
