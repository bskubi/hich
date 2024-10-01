include {QCReads} from './qcHicReads.nf'
include {emptyOnLastStep; pack} from './extraops.nf'

process HichFragtag {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null},
               mode: params.general.publish.mode
    container "bskubi/hich:latest"
    label 'doJobArray'
    label 'pairs'

    input:
    tuple val(id), path(pairs), path(fragmentIndex)

    output:
    tuple val(id), path("${id}_fragtag.pairs.gz")

    shell:
    "hich fragtag ${fragmentIndex} ${id}_fragtag.pairs.gz ${pairs}"

    stub:
    "touch ${id}_fragtag.pairs.gz"
}

def hasFragmentIndexName = {it.get("fragmentIndex").toString().trim().length() > 0}

def fragmentIndexExists = {hasFragmentIndexName(it) && file(it.fragmentIndex).exists()}

workflow TagRestrictionFragments {
    take:
        samples

    main:

    samples
        | filter{fragmentIndexExists(it) && (it.pairs || it.latestPairs)}
        | map{tuple(it.id, it.pairs, it.fragmentIndex)}
        | HichFragtag
        | map{[id:it[0], fragPairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{result}
    pack(samples, result) | set{samples}

    if ("TagFragments" in params.general.get("qcAfter")) {
        samples = QCReads(samples, "TagFragments")
    }

    samples = emptyOnLastStep("TagFragments", samples)

    emit:
        samples
}