include {pack; skip; coalesce; rows; columns; groupHashMap; isTechrep; isBiorep; isCondition; constructIdentifier; emptyOnLastStep; aggregateLevelLabel} from './extraops.nf'
include {Setup} from './setup.nf'

/*
    Aggregation can involve downsampling, deduplicating, and merging samples, all
    steps being optional. In the config, the params.aggregate section defines
    aggregateProfiles that define how samples should be downsampled, deduplicated and merged
    at each level. An aggregateLevel is techreps, bioreps, or conditions.

    Overall workflow
        Aggregate[Level] -> CreateAggregateSampleProfiles -> Downsample -> Deduplicate -> Merge to [Level + 1]

    Individual downsampling
    -----------------------
    Downsampling can be simply to a specific fraction of original reads or exact number of reads.

    Group downsampling
    ------------------
    Downsampling can also partition the reads in each sample according to conjuncts over PairsSegment
    attributes specified by the user, such as chrom1, chrom2, and pair_type, as well as user-defined
    cis-strata. The fraction of reads in each read partition block (outcome in read sample space)
    can be computed, averaged (over samples not tagged as "outliers") and then used as a target
    distribution for the input samples, which will be minimally downsampled to conform to this
    target mean non-outlier distribution and then further downsampled to the size of the smallest
    member of the comparison group. This exploits a triplicate replicate structure to filter
    out technical noise while controlling for sample-sample coverage differences.

    Group downsampling operates only within an aggregationProfile.

    Downsampling implementation
    ---------------------------
    Downsampling begins by taking a list of conjuncts and calling hich stats, which counts the number
    of reads for each outcome (conjunct combination).

    hich stats-aggregate is then called on groups of stats files (which can be groups of size 1 if no real
    group downsampling is scheduled for a particular aggregation profile). This handles computing the target
    number of per-sample reads to retain and saves the results as a new batch of stats files, one for each
    sample.

    Often, reservoire sampling is used to sample a specific number of events in one pass through a dataset, but
    since this workflow intrinsically requires counting the starting number of reads, we can instead use a much
    simpler and lower-memory algorithm published by Donald Knuth in The Art of Computer Programming, Vol. 2, "Algorithm S".
    The original and target distribution of read counts is passed for each individual sample to hich downsample,
    which samples the specific number of reads and outputs to a new pairs file.

    Level-specific downsampling
    ---------------------------
    Each aggregation profile can set the approach to downsampling differently for each aggregateLevel (techrep, biorep, condition)

    Deduplication
    -------------
    Options include
        No deduplication
        Computing wobble as either
            W = sum(abs(r1.pos1 - r2.pos1) + abs(r1.pos2 - r2.pos2))
            W = max(abs(r1.pos1 - r2.pos1) + abs(r1.pos2 - r2.pos2))
            Then requiring W < T for arbitrary non-negative T to call two reads as duplicates
        Backend selection
            KD Tree backend (default for non-single cell) is transitive, meaning that if
            r1 and r2 are duplicates and r2 and r3 are duplicates, then r1 and r3 are duplicates.

            Cython backend (default for single-cell for performance reasons) is non-transitive,
            meaning that in the above scenario r1 and r3 are not necessarily duplicates
        Single-cell deduplication
            If a sample is marked "isSingleCell", it's expected that as it's parsed to pairs format,
            it contains a cellID column containing a barcode corresponding to the cell from which it originates.
            This will be used as an extra column that must match during deduplication for two reads to be
            marked as duplicates.
    Samples are individually deduplicated only after they have (optionally) been merged, meaning that there will
    be no redundant deduplication. If T1 and T2 are merged to B1, then T1 and T2 will only subsequently be deduplicated.
    Techreps are deduplicated on the AggregateTechreps step, bioreps on the AggregateBioreps step, conditions on
    AggregateConditions, even though the bioreps and conditions are produced via a merge on the previous step.

    Merge
    -----
    Merges operates only within an aggregationProfile and yield a sample promoted one aggregateLevel from the inputs.
    It takes place on the lower level (i.e. for techreps -> bioreps, during AggregateTechreps).
*/


/*
    Compute original read counts over conjuncts defined for the current agg level and profile with hich stats.
*/
process HichStats {
    container "bskubi/hich:latest"
    label 'pairs'

    input:
    tuple val(id), path(pairs), val(conjuncts), val(cisStrata)

    output:
    tuple val(id), path(stats)

    shell:
    stats = "${id}.from.stats.tsv"
    conjunctsArg = conjuncts ? "--conjuncts '${conjuncts.join(' ')}'" : ""
    cisStrataArg = cisStrata ? "--cis-strata '${cisStrata.join(' ')}'" : ""
    "hich stats ${conjunctsArg} ${cisStrataArg} --output ${stats} ${pairs}"

    stub:
    stats = "${id}.stats.tsv"
    "touch ${stats}"
}

/*
    Compute target read counts over groups and their common conjuncts defined for the current agg level and profile
    with hich stats-aggregate
*/
process HichStatsAggregate {
    //container "hich:latest"
    //container "bskubi/hich:latest"
    label 'pairs'

    input:
    tuple val(ids), path(stats), val(outliers)

    output:
    tuple val(ids), path(targetStats)

    shell:
    targetStats = stats.collect{"aggregate_${it}"}
    outlier_stats = stats.findAll { stats.indexOf(it) in outliers.findIndexValues { it } }
    outlier_params = outlier_stats.collect{"--outlier ${it}"}.join(" ")
    "hich stats-aggregate --to-group-mean --to-group-min --prefix aggregate_ ${outlier_params} ${stats.join(' ')}"
    
    stub:
    targetStats = stats.collect{"aggregate_${it}"}
    "touch ${targetStats.join(' ')}"
}

/*
    Now that the original and target distribution is computed for each sample, individually downsample them to the target
    distribution with hich downsample.
*/
process HichDownsample {
    //container "bskubi/hich:latest"
    label 'pairs'

    input:
    tuple val(id), path(fullPairs), path(statsFrom), path(statsTo), val(conjuncts), val(cisStrata), val(toSize)

    output:
    tuple val(id), path(downsampledPairs)

    shell:
    downsampledPairs = "${id}.downsampled.pairs.gz"
    conjunctsArg = conjuncts ? "--conjuncts '${conjuncts.join(' ')}'" : ""
    cisStrataArg = cisStrata ? "--cis-strata '${cisStrata.join(' ')}'" : ""
    toSizeArg = toSize ? "--to-size ${toSize}" : ""

    "hich downsample ${conjunctsArg} ${cisStrataArg} --orig-stats ${statsFrom} --target-stats ${statsTo} ${toSizeArg} ${fullPairs} ${downsampledPairs}"

    stub:
    downsampledPairs = "${id}.downsampled.pairs.tsv"
    "touch ${downsampledPairs}"
}

/*
    Deduplicate pairs files
*/

process PairtoolsDedup {
    publishDir params.general.publish.parse ? params.general.publish.dedup : "results",
               saveAs: {params.general.publish.dedup ? it : null},
               mode: params.general.publish.mode
    conda "bioconda::pairtools bioconda::samtools"
    container "bskubi/hich:latest"
    label 'pairs'
    cpus 8

    input:
    tuple val(id), path(pairs), val(singleCell), val(maxMismatch), val(method), val(pairtoolsDedupParams)

    output:
    tuple val(id), path(deduplicated)

    shell:
    deduplicated = "${id}.dedup.pairs.gz"
    
    if (!pairtoolsDedupParams) pairtoolsDedupParams = []
    if (singleCell) pairtoolsDedupParams += ["--extra-col-pair cellID cellID --backend cython --chunksize 1000000"]
    if (maxMismatch != null) pairtoolsDedupParams += ["--max-mismatch ${maxMismatch}"]
    if (method != null) pairtoolsDedupParams += ["--method ${method}"]

    cmd = "pairtools dedup --output ${deduplicated} ${pairtoolsDedupParams.join(' ')}  --nproc-in ${task.cpus} --nproc-out ${task.cpus} ${pairs}"
    cmd

    stub:
    deduplicated = "${id}.dedup.pairs.gz"
    "touch ${deduplicated}"
}

/*
    Merge pairs files to one level higher (techreps -> bioreps -> conditions)
*/
process PairtoolsMerge {
    container "bskubi/hich:latest"
    label 'pairs'
    cpus 8
    
    input:
    tuple val(id), path(to_merge)

    output:
    tuple val(id), path(merged)

    shell:
    merged = "${id}.merged.pairs.gz"
    "pairtools merge --output ${merged}  --nproc-in ${task.cpus} --nproc-out ${task.cpus} ${to_merge.join(' ')}"

    stub:
    merged = "${id}.merged.pairs.gz"
    "touch ${merged}"
}


/*
    This process clones raw samples at the current agg level and assigns them
    to an agg profile, applying the params for that profile to the sample.

    The nextflow.config file may have a params.aggregate {} section defining
    aggregate profiles. If so, create duplicates for each sample, one for each
    profile, and label the sample with the aggregate profile, keeping it separate
    from samples for other agg profiles going forward.

    We generalize the main aggregation workflow steps but have an initial level-specific
    workflow that assigns level-specific params with general names in a levelParams map
    of [generalKey: levelSpecificKey].

    levelSamples -- the samples filtered for the current level (i.e. techreps, bioreps, conditions)
    levelParams -- the map from general aggregation keywords to the level-specific keywords contained in the sample map

    Note that this calls no processes.
*/
workflow CreateAggregateSampleProfiles {
    take:
    levelSamples
    levelParams

    main:

    // Store current samples as the result in case there are no aggregation profiles defined
    levelSamples | set{result}

    if (params.containsKey("aggregate")) {
        /* The user may input samples from multiple aggregation levels. Raw samples not assigned
            an aggregation profile need to be cloned for each aggregation level so this selects them
            for this cloning process.
        */
        levelSamples | branch{
            YES: it.aggregateProfileName == null
            NO: true
        } | set{raw}

        /*
            From the RAW SAMPLES ONLY, create new samples associated with each aggregation profile

            First, reshape the profiles to a columnar format like [profileName: [list of profile names], profileParams: [list of profile param maps]]
            Then convert to a row format and integrate with the samples, ultimately resulting in samples with new parameters:
            [original sample attributes] + [aggregateProfileName: profileName, profileParams: profileParams]

            We give each newly created sample profile-specific id.
        */
        aggregateProfiles = channel.of(params.aggregate)
        | map{[profileName: it.keySet().toList(), profileParams: it.values().toList()]}
        | map{rows(it)}
        | flatten
        | combine(raw.YES)
        | map {
            profile = it[0]
            sample = it[1]
            sample += [aggregateProfileName: profile.profileName] + profile.profileParams
            sample
        }
        | map{
            downsampleToMeanDistributionGroup = it.subMap((levelParams.downsampleToMeanDistribution)) ?: null
            it + [(levelParams.downsampleToMeanDistribution): downsampleToMeanDistributionGroup] + [id: constructIdentifier(it)]
        }
        | set{result}

        /*
            Typically, we will expect the user to drop the raw samples, specifying an explicit raw-like
            aggregation profile if they want that. However, the user can also set a keepUnaggregated
            flag in order to retain the raw samples.
        */

        if (params.containsKey("keepUnaggregated")) {
            result | concat(raw.NO) | set{result}
        }
    }

    emit:
    result
    levelParams
}

workflow DownsamplePairs {
    take:
    levelSamples
    levelParams

    main:
    // levelParams gets converted to a DataflowVariable when emitted from the previous workflow,
    // so extract the LinkedHashMap value so we can use it more conveniently
    levelParams = levelParams.value

    levelSamples | branch {
        YES: it.aggregateProfileName != null && (it.get(levelParams.downsampleToMeanDistribution) != null || !(it.get(levelParams.downsampleToSize) in [1.0, null]))
        NO: true
    } | set{downsample}

    // Compute from stats
    downsample.YES
        | map{tuple(it.id, it.latestPairs, it.get(levelParams.readConjuncts), it.get(levelParams.cisStrata))}
        | HichStats
        | map{[id:it[0], (levelParams.downsampleStatsFrom):it[1]]}
        | set{downsampleStatsFromResult}
    pack(downsample.YES, downsampleStatsFromResult) | set {fromStatsCalculated}

    // Compute to stats
    groupHashMap(fromStatsCalculated, [levelParams.downsampleToMeanDistribution, 'aggregateProfileName'])
        | map{columns(it, ["dropNull":true])}
        | map{tuple(it.id, it.get(levelParams.downsampleStatsFrom), it.outlier)}
        | HichStatsAggregate
        | map{[id:it[0], (levelParams.downsampleStatsTo):it[1]]}
        | map{rows(it)}
        | flatten
        | set {downsampleToGroupMinResult}
    pack(fromStatsCalculated, downsampleToGroupMinResult) | set{toGroupMinResult}

    // Downsample
    toGroupMinResult
        | map{
            tuple(it.id, it.latestPairs, it.get(levelParams.downsampleStatsFrom), it.get(levelParams.downsampleStatsTo),
              it.get(levelParams.readConjuncts), it.get(levelParams.cisStrata), it.get(levelParams.downsampleToSize))}
        | HichDownsample
        | map{[id:it[0], (levelParams.downsamplePairs):it[1], latest:it[1], latestPairs:it[1]]}
        | set{downsampledResult}
    pack(toGroupMinResult, downsampledResult) | set {downsampled}

    downsampled | concat(downsample.NO) | set{result}

    emit:
    result
    levelParams
}

workflow MergePairs {
    take:
    levelSamples
    levelParams

    main:

    // levelParams gets converted to a DataflowVariable when emitted from the previous workflow,
    // so extract the LinkedHashMap value so we can use it more conveniently
    levelParams = levelParams.value

    levelSamples | branch {
        YES: it.includeInMerge != false && it.get(levelParams.doMerge)
        NO: true}
    | set{merge}

    // Merge
    groupHashMap(merge.YES, levelParams.mergeGroupIdentifiers)
        | map{coalesce(columns(it, ["dropNull":true]))}
        | map{tuple(constructIdentifier(coalesce(it, "_drop")), it.latestPairs)}
        | PairtoolsMerge
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{fromMerge}

    groupHashMap(merge.YES, levelParams.mergeGroupIdentifiers)
        | map{coalesce(columns(it, ["dropNull":true]), '_drop')}
        | map{it += [id:constructIdentifier(it)]}
        | set{inheritedMergeAttributes}
    pack(inheritedMergeAttributes, fromMerge) | set{merged}

    merged | map {it += ["aggregateLevel" : aggregateLevelLabel(it)]} | concat(levelSamples) | set{result}

    emit:
    result
    levelParams
}

workflow DeduplicatePairs {
    take:
    levelSamples
    levelParams

    main:
    
    // levelParams gets converted to a DataflowVariable when emitted from the previous workflow,
    // so extract the LinkedHashMap value so we can use it more conveniently
    levelParams = levelParams.value

    levelSamples
        | branch {
            YES: !it.alreadyDeduplicated && levelParams.levelFilter(it) && it.get(levelParams.doDedup)
            NO: true
    } | set{deduplicate}

    deduplicate.YES
        | map{tuple(it.id, it.latestPairs, it.dedupSingleCell, it.dedupMaxMismatch, it.dedupMethod, it.pairtoolsDedupParams)}
        | PairtoolsDedup
        | map{[id:it[0], dedupPairs:it[1], latest:it[1], latestPairs:it[1], alreadyDeduplicated:true]}
        | set{dedupResult}
    pack(deduplicate.YES, dedupResult) | concat(deduplicate.NO) | set{result}

    emit:
    result
}

workflow AggregateTechreps {
    take:
    samples

    main:

    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'max',
        cisStrata: 'techrepCisStrata',
        readConjuncts: 'techrepReadConjuncts',
        downsampleToMeanDistribution: 'techrepDownsampleToMeanDistribution',
        downsampleToSize: 'techrepDownsampleToSize',
        downsampleStatsFrom: 'techrepDownsampleStatsFrom',
        downsampleStatsTo: 'techrepDownsampleStatsTo',
        downsamplePairs: 'techrepDownsamplePairs',
        doMerge: 'mergeTechrepToBiorep',
        doDedup: 'techrepDedup',
        mergeGroupIdentifiers: ['condition', 'biorep', 'aggregateProfileName'],
        levelFilter: {it.condition && it.biorep && it.techrep}
    ]

    samples | branch {
        techrep: !skip("aggregate") && !skip("aggregateTechreps") && isTechrep(it) && (it.pairs || it.latestPairs)
        other: true
    } | set{sampleType}

    CreateAggregateSampleProfiles(sampleType.techrep, aggregateParams)
    | DownsamplePairs
    | MergePairs
    | DeduplicatePairs
    | concat(sampleType.other)
    | set{samples}

    samples = emptyOnLastStep("aggregateTechreps", samples)

    emit:
    samples
}

workflow AggregateBioreps {
    take:
    samples

    main:

    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'max',
        cisStrata: 'biorepCisStrata',
        readConjuncts: 'biorepReadConjuncts',
        downsampleToMeanDistribution: 'biorepDownsampleToMeanDistribution',
        downsampleToSize: 'biorepDownsampleToSize',
        downsampleStatsFrom: 'biorepDownsampleStatsFrom',
        downsampleStatsTo: 'biorepDownsampleStatsTo',
        downsamplePairs: 'biorepDownsamplePairs',
        doMerge: 'mergeBiorepToCondition',
        doDedup: 'biorepDedup',
        mergeGroupIdentifiers: ['condition', 'aggregateProfileName'],
        levelFilter: {it.condition && it.biorep && !it.techrep}
    ]

    samples | branch {
        biorep: !skip("aggregate") && !skip("aggregateBioreps") && isBiorep(it) && (it.pairs || it.latestPairs)
        other: true
    } | set{sampleType}

    CreateAggregateSampleProfiles(sampleType.biorep, aggregateParams)
    | DownsamplePairs
    | MergePairs
    | DeduplicatePairs
    | concat(sampleType.other)
    | set{samples}

    samples = emptyOnLastStep("aggregateBioreps", samples)

    emit:
    samples
}

workflow AggregateConditions {
    take:
    samples

    main:

    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'max',
        cisStrata: 'conditionCisStrata',
        readConjuncts: 'conditionReadConjuncts',
        downsampleToMeanDistribution: 'conditionDownsampleToMeanDistribution',
        downsampleToSize: 'conditionDownsampleToSize',
        downsampleStatsFrom: 'conditionDownsampleStatsFrom',
        downsampleStatsTo: 'conditionDownsampleStatsTo',
        downsamplePairs: 'conditionDownsamplePairs',
        doMerge: 'mergeCondition',
        doDedup: 'conditionDedup',
        mergeGroupIdentifiers: ['aggregateProfileName'],
        levelFilter: {it.condition && !it.biorep && !it.techrep}
    ]

    samples | branch {
        condition: !skip("aggregate") && !skip("aggregateConditions") && isCondition(it) && (it.pairs || it.latestPairs)
        other: true
    } | set{sampleType}

    CreateAggregateSampleProfiles(sampleType.condition, aggregateParams)
    | DownsamplePairs
    | MergePairs
    | DeduplicatePairs
    | concat(sampleType.other)
    | set{samples}

    samples = emptyOnLastStep("aggregateConditions", samples)

    emit:
    samples
}
