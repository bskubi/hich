include {QCReads} from './qcHicReads.nf'
include {AssignParams} from './assignParams.nf'
include {transpack} from './extraops.nf'

process PairtoolsMerge {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null},
               mode: params.general.publish.mode
    conda "pairtools"
    container "bskubi/hich:latest"

    input:
    tuple val(id), path(samples)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    samples = (samples.getClass() == nextflow.processor.TaskPath) ? samples : samples.join(" ")
    "pairtools merge --output ${id}.pairs.gz ${samples}"

    stub:
    "touch ${id}.pairs.gz"
}

def label(hashmap, label) {
    hashmap.containsKey(label) &&
    hashmap.get(label) != null &&
    hashmap.get(label).toString().trim().length() >= 1
}

def isTechrep = {map -> label(map, "techrep") && label(map, "biorep") && label(map, "condition")}
def isBiorep = {map -> !label(map, "techrep") && label(map, "biorep") && label(map, "condition")}
def isCondition = {map -> !label(map, "techrep") && !label(map, "biorep") && label(map, "condition")}

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
    if (isTechrep(hashmaps[0])) {
        merge.remove("techrep")
        merge.id = "${merge.condition}_${merge.biorep}"
    } else if (isBiorep(hashmaps[0])) {
        merge.remove("biorep")
        merge.id = merge.condition
    } else {
        error "merge ${merge.id} is not a condition or biorep:\n${merge}"
    }
    return merge
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