include {JoinProcessResults} from './joinProcessResults.nf'
include {QCReads} from './qcHicReads.nf'
include {AssignParams} from './assignParams.nf'

process Merge {
    publishDir params.general.publish.fragtag ? params.general.publish.fragtag : "results",
               saveAs: {params.general.publish.fragtag ? it : null}
            
    container "bskubi/pairtools:1.1.0"

    input:
    tuple val(id), path(samples)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    samples = (samples.getClass() == nextflow.processor.TaskPath) ? samples : samples.join(" ")
    "pairtools --version > version.txt && pairtools merge --output ${id}.pairs.gz ${samples}"
}

workflow TechrepsToBioreps {
    take:
        samples

            main:

        def isTechrep = {it.get("is_techrep", "").toString().trim() in ["", "true", true, 1]}
        def hasStructure = {it.get("condition").length() > 0
                            && it.get("biorep").length() > 0}
        def ensureStructure = {
            if (isTechrep(it)) {it.is_techrep = true}
            if (isTechrep(it) && !hasStructure(it)) {
                it.biorep = it.sample_id
                it.condition = it.sample_id
            }
            it
        }

        samples
            | filter{isTechrep(it)}
            | map{ensureStructure(it)}
            | map{tuple(it.subMap("condition", "biorep"), it)}
            | groupTuple
            | map {
                biorep, techreps ->
                biorep += techreps[0].findAll{entry -> techreps.every {map -> map[entry.key] == entry.value}}
                biorep.latest = techreps.collect{it.latest}
                biorep.is_biorep = true
                biorep.id = "${biorep.condition}_${biorep.biorep}"
                biorep
            }
            | AssignParams
            | set{to_merge}

        to_merge = JoinProcessResults(
            Merge,
            [to_merge],
            ["id", "latest"],
            ["id", "biorep_merge_pairs"],
            ["id"],
            false,
            "biorep_merge_pairs"
        )

        to_merge
            | concat(samples)
            | set {samples}

        if ("TechrepsToBioreps" in params.general.get("qc_after")) {
            QCReads(samples, "TechrepsToBioreps")
        }

    emit:
        samples
}

workflow BiorepsToConditions {
    take:
        samples

    main:

        def isBiorep = {it.containsKey("is_biorep") && it.is_biorep}
        def hasStructure = {it.get("condition").length() > 0
                            && it.get("biorep").length() > 0}
        def ensureStructure = {
            if (isBiorep(it)) {it.is_biorep = true}
            if (isBiorep(it) && !hasStructure(it)) {
                it.biorep = it.sample_id
                it.condition = it.sample_id
            }
            it
        }

        samples
            | filter{isBiorep(it)}
            | map{ensureStructure(it)}
            | map{tuple(it.subMap("condition"), it)}
            | groupTuple
            | map {
                condition, bioreps ->
                condition += bioreps[0].findAll{entry -> bioreps.every {map -> map[entry.key] == entry.value}}
                condition.latest = bioreps.collect{it.latest}
                condition.id = "${condition.condition}"
                condition.is_condition = true
                condition
            }
            | AssignParams
            | set{to_merge}


        to_merge = JoinProcessResults(
            Merge,
            [to_merge],
            ["id", "latest"],
            ["id", "condition_merge_pairs"],
            ["id"],
            false,
            "condition_merge_pairs"
        )

        to_merge
            | concat(samples)
            | set {samples}

        if ("BiorepsToConditions" in params.general.get("qc_after")) {
            QCReads(samples, "BiorepsToConditions")
        }

    emit:
        samples
}