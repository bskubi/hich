include {QCReads} from './qcHicReads.nf'
include {Setup} from './setup.nf'
include {transpack; isTechrep; isBiorep; isCondition; emptyOnLastStep} from './extraops.nf'

process PairtoolsMerge {
    publishDir params.general.publish.merge ? params.general.publish.merge : "results",
               saveAs: {params.general.publish.merge ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::pairtools"
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

// Get a hashmap containing only key:value pairs identical across all HashMaps.
def getSameValueSubset = {
    hashmaps ->

    // Implementation:
    // Start with the first HashMap and iterate through each of its keys
    // For each key:
    //     Keep only if its value in the first HashMap is the same for every HashMap
    // Return the subset of the hashmap for the kept keys only
    hashmaps[0].findAll {
        item ->
        
        hashmaps.every {
            map -> map[item.key] == item.value
        }
    }
}

// Drop the non-shared subset of parameters across all HashMaps
// Except latest, which are combined to a list.
// Remove the "techrep" or "biorep" key as appropriate.
// Then create a composite id.
def mergeGroupParams = {
    hashmaps ->

    // Keep key:value pairs that are identical across all hashmaps.
    // Drop the rest of the key:value pairs.
    def merge = getSameValueSubset(hashmaps)

    // Create a list of the "latest" pairs for every hashmap
    merge.latest = []
    hashmaps.each {map -> merge.latest.add(map.get("latest"))}
    merge.latest = merge.latest.clone().sort()

    // Edit the indicators for being a techrep/biorep as necessary and
    // create a composite id based on the merged sampleType
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
        case "TechrepsToBioreps":
            groups = ["condition", "biorep"]
            result = "biorep_merge_pairs"
            break
        case "BiorepsToConditions":
            groups = ["condition"]
            result = "condition_merge_pairs"
            break
        default:
            error "Invalid merge sampleType '${sampleType}'"
    }

    // Filter for the samples to merge
    // Group them according to their similars (i.e. same condition/biorep)
    // Then drop the non-shared keys, put the "latest" values into a list,
    // drop the "techrep" or "biorep" keys as appropriate, and create a composite ID
    // Finally, apply the Setup parameters based on the composite ID
    samples
        | filter{
            switch(sampleType) {
                case "TechrepsToBioreps":
                    isTechrep(it)
                    break
                case "BiorepsToConditions":
                    isBiorep(it)
                    break
                default:
                    error "Invalid merge sampleType '${sampleType}'"
            }
        }
        | map{tuple(it.subMap(*groups), it)}
        | groupTuple
        | map {group, samples -> mergeGroupParams(samples)}
        | Setup
        | set{to_merge}

    // The "latest" key now contains a list of the pairs to merge
    // Store the result in biorep_merge_pairs or condition_merge_pairs
    // as well as in "latest"
    to_merge = transpack(
        PairtoolsMerge,
        [to_merge],
        ["id", "latest"],
        ["id", result],
        ["latest":result]
    )
    
    // Concatenate the merged samples to the rest of the samples
    to_merge
        | concat(samples)
        | set {samples}

    if (sampleType in params.general.get("qcAfter")) {
        QCReads(samples, sampleType)
    }

    emit:
    samples
}

workflow TechrepsToBioreps {
    take:
    samples

    main:

    samples = Merge(samples, "TechrepsToBioreps")

    samples = emptyOnLastStep("TechrepsToBioreps", samples)

    emit:
    samples
}

workflow BiorepsToConditions {
    take:
    samples

    main:
    samples = Merge(samples, "BiorepsToConditions")
    samples = emptyOnLastStep("BiorepsToConditions", samples)

    emit:
    samples
}