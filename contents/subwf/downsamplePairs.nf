include {transpack; isTechrep; isBiorep; isCondition; emptyOnLastStep} from './extraops.nf'

process HichDownsample {
    input:
    tuple val(id), path(pairs), path(orig_stats), path(target_stats)

    output:
    tuple val(id), path("${id}_downsampled.pairs.gz")

    shell:
    "hich downsample --orig_stats ${orig_stats} --target_stats ${target_stats} ${pairs} ${id}_downsampled.pairs.gz"

    stub:
    "touch ${id}_downsampled.pairs.gz"
}

process HichStats {
    input:
    tuple val(id), path(pairs), val(hich_stats_params)

    output:
    tuple val(id), path("${id}.stats.tsv")

    shell:
    hich_stats_params = hich_stats_params.join(" ")
    "hich stats ${hich_stats_params} ${pairs} ${id}.stats.tsv"
}

workflow Downsample {
    take:
    samples
    sampleType

    main:
    // Get the sample traits to group by ("groups") and the sample trait name
    // to save the resulting pairs file as ("downsample_techreps_pairs")
    switch(sampleType) {
        case "DownsampleTechreps":
            groups = ["condition", "biorep"]
            pairs = "downsample_techreps_pairs"
            orig_stats = "techreps_stats"
            target_stats = "downsample_techreps_target"
            break
        case "DownsampleBioreps":
            groups = ["condition"]
            pairs = "downsample_bioreps_pairs"
            orig_stats = "bioreps_stats"
            target_stats = "downsample_bioreps_target"
            break
        case "DownsampleConditions":
            groups = []
            pairs = "downsample_conditions_pairs"
            orig_stats = "conditions_stats"
            target_stats = "downsample_techreps_target"
            break
        default:
            error "Invalid merge sampleType '${sampleType}'"
    }

    // First, filter the samples for only techreps, bioreps, or conditions as appropriate.
    samples
        | filter{
            switch(sampleType) {
                case "DownsampleTechreps":
                    isTechrep(it)
                    break
                case "DownsampleBioreps":
                    isBiorep(it)
                    break
                case "DownsampleConditions":
                    isCondition(it)
                    break
                default:
                    error "Invalid merge sampleType '${sampleType}'"
            }
        }
        | map{it[orig_stats]}
        | set{to_downsample}

    // Second, call stats on to_downsample with transpack
    to_downsample = transpack(
        HichStats,
        [to_downsample],
        ["id", "latest", "hichStatsParams"],
        ["id", orig_stats],
        [],
        "id"
    )

    // Then put the traits ("condition", "stats", "none") as the first item in a tuple with the items themselves as the second
    // Then do a groupTuple on the first item to group the filtered samples with their similars
    // (other techreps in the same condition + biorep). This will have the form [[grouping], [list of HashMap items]]
    to_downsample
        | map{tuple(it.subMap(*groups), it)}
        | groupTuple
    // Filter for appropriate sampleType
    // Call stats and store as original distribution
    // Call stats-aggregate and store as target distribution
    // Call downsample and store as latest + "downsampled"

    emit:
    samples
}

workflow DownsampleConditions {
    take:
    samples

    main:
    samples = Downsample(samples, "DownsampleConditions")

    samples = emptyOnLastStep("DownsampleConditions", samples)

    emit:
    samples
}


workflow DownsampleBioreps {
    take:
    samples

    main:
    samples = Downsample(samples, "DownsampleBioreps")

    samples = emptyOnLastStep("DownsampleBioreps", samples)


    emit:
    samples
}

workflow DownsampleTechreps {
    take:
    samples

    main:
    samples = Downsample(samples, "DownsampleTechreps")

    samples = emptyOnLastStep("DownsampleTechreps", samples)
    
    emit:
    samples
}