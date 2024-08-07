include {QCReads} from './qcHicReads.nf'
include {AssignParams} from './assignParams.nf'
include {transpack; isTechrep; isBiorep; isCondition} from './extraops.nf'

process PairtoolsMerge {
    publishDir params.general.publish.merge ? params.general.publish.merge : "results",
               saveAs: {params.general.publish.merge ? it : null},
               mode: params.general.publish.mode
    conda "pairtools"
    container "bskubi/hich:latest"

    input:
    tuple val(id), path(samples)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    def to_merge = (samples.getClass() == nextflow.processor.TaskPath) ? samples : samples.join(" ")
    "pairtools merge --output ${id}.pairs.gz ${to_merge}"

    stub:
    "touch ${id}.pairs.gz"
}

def mergeHashmaps = {
    hashmaps ->

    hashmaps[0].findAll {
        item ->
        
        hashmaps.every {
            map -> map[item.key] == item.value
        }
    }
}

def mergeGroupParams = {
    hashmaps ->
    def merge = mergeHashmaps(hashmaps)
    merge.latest = []
    hashmaps.each {map -> merge.latest.add(map.get("latest"))}
    merge.latest = merge.latest.clone().sort()
    if (isTechrep(hashmaps[0])) {
        merge.remove("techrep")
        merge.id = "${merge.condition}_${merge.biorep}"
    } else if (isBiorep(hashmaps[0])) {
        merge.remove("biorep")
        merge.id = merge.condition
    } else {
        error "merge ${merge.id} is not a condition or biorep:\n${merge}"
    }
    merge
}

workflow Merge {
    take:
    samples
    sampleType

    main:

    switch(sampleType) {
        case "TechrepsToBiorep":
            filterFunction = isTechrep
            groups = ["condition", "biorep"]
            result = "biorep_merge_pairs"
            break
        case "BiorepsToCondition":
            filterFunction = isBiorep
            groups = ["condition"]
            result = "condition_merge_pairs"
            break
        default:
            error "Invalid merge sampleType '${sampleType}'"
    }

    samples
        | filter{filterFunction.call(it)}
        | map{tuple(it.subMap(*groups), it)}
        | groupTuple
        | map {group, samples -> mergeGroupParams(samples)}
        | AssignParams
        | set{to_merge}

    to_merge = transpack(
        PairtoolsMerge,
        [to_merge],
        ["id", "latest"],
        ["id", result],
        ["latest":result]
    )
    
    to_merge
        | concat(samples)
        | set {samples}

    if (sampleType in params.general.get("qc_after")) {
        QCReads(samples, sampleType)
    }

    emit:
    samples
}

workflow TechrepsToBioreps {
    take:
    samples

    main:

    samples = Merge(samples, "TechrepsToBiorep")

    emit:
    samples
}

workflow BiorepsToConditions {
    take:
    samples

    main:
    samples = Merge(samples, "BiorepsToCondition")

    emit:
    samples
}