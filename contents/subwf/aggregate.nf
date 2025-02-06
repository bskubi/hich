include {withLog; stubLog; pack; skip; coalesce; rows; columns; groupHashMap; isTechrep; isBiorep; isCondition; constructIdentifier; emptyOnLastStep; aggregateLevelLabel} from './extraops.nf'
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
    container params.general.hichContainer
    label 'pairs'

    input:
    tuple val(id), path(pairs), val(conjuncts), val(cisStrata)

    output:
    tuple val(id), path(stats)

    shell:
    // Create name of stats file
    stats = "${id}.from.stats.tsv"

    // Format options --conjuncts and --cis-strata
    conjunctsArg = conjuncts ? "--conjuncts '${conjuncts.join(',')}'" : ""
    cisStrataArg = cisStrata ? "--cis-strata '${cisStrata.join(',')}'" : ""
    
    // Format the hich stats command
    cmd = "hich stats ${conjunctsArg} ${cisStrataArg} --output '${stats}' '${pairs}'"

    // Store the parameters and command for this call to the HichStats process
    logMap = [task: "HichStats", input: [id: id, conjuncts: conjunctsArg, cisStrata: cisStrataArg, pairs: pairs], output: stats]

    // Log the command and execute
    withLog(cmd, logMap)

    stub:
    stats = "${id}.from.stats.tsv"
    stub = "touch ${stats}"
    conjunctsArg = conjuncts ? "--conjuncts '${conjuncts.join(',')}'" : ""
    cisStrataArg = cisStrata ? "--cis-strata '${cisStrata.join(',')}'" : ""
    cmd = "hich stats ${conjunctsArg} ${cisStrataArg} --output '${stats}' '${pairs}'"
    logMap = [task: "HichStats", input: [id: id, conjuncts: conjunctsArg, cisStrata: cisStrataArg, pairs: pairs], output: stats]
    stubLog(stub, cmd, logMap)
}

/*
    Compute target read counts over groups and their common conjuncts defined for the current agg level and profile
    with hich stats-aggregate
*/
process HichStatsAggregate {
    container params.general.hichContainer
    label 'pairs'

    input:
    tuple val(ids), path(stats), val(outliers)

    output:
    tuple val(ids), path(targetStats)

    shell:
    
    if (stats == null || !(stats.metaClass.respondsTo(stats, "indexOf"))) {stats = [stats]}
    stats = stats.collect { "${it}" }

    if (outliers == null || !(outliers.metaClass.respondsTo(outliers, "findIndexValues"))) {outliers = [outliers]}
    outliers = outliers.collect { "'${it}'" }

    targetStats = stats.collect{"aggregate_${it}"}
    outlier_stats = stats.findAll { stats.indexOf(it) in outliers.findIndexValues { it } }
    outlier_params = outlier_stats.collect{"--outlier ${it}"} ?: ""
    cmd = "hich stats-aggregate --to-group-mean --to-group-min --prefix aggregate_ ${outlier_params} ${stats.join(' ')}"
    logMap = [task: "HichStatsAggregate", input: [id: ids, stats: stats, outliers: outliers], output: targetStats]
    withLog(cmd, logMap)

    stub:
    stub = "touch ${targetStats.join(' ')}"
    if (stats == null || !(stats.metaClass.respondsTo(stats, "indexOf"))) {stats = [stats]}
    stats = stats.collect { "${it}" }

    if (outliers == null || !(outliers.metaClass.respondsTo(outliers, "findIndexValues"))) {outliers = [outliers]}
    outliers = outliers.collect { "'${it}'" }

    targetStats = stats.collect{"aggregate_${it}"}
    outlier_stats = stats.findAll { stats.indexOf(it) in outliers.findIndexValues { it } }
    outlier_params = outlier_stats.collect{"--outlier ${it}"} ?: ""
    cmd = "hich stats-aggregate --to-group-mean --to-group-min --prefix aggregate_ ${outlier_params} ${stats.join(' ')}"
    logMap = [task: "HichStatsAggregate", input: [id: ids, stats: stats, outliers: outliers], output: targetStats]
    stubLog(stub, cmd, logMap)

}

/*
    Now that the original and target distribution is computed for each sample, individually downsample them to the target
    distribution with hich downsample.
*/
process HichDownsample {
    container params.general.hichContainer
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
    cmd = "hich downsample ${conjunctsArg} ${cisStrataArg} --orig-stats ${statsFrom} --target-stats ${statsTo} ${toSizeArg} ${fullPairs} ${downsampledPairs}"
    logMap = [
        task: "HichDownsample",
        output: downsampledPairs,
        input: [
            id: id,
            fullPairs: fullPairs,
            statsFrom: statsFrom,
            statsTo: statsTo,
            conjuncts: conjuncts,
            cisStrata: cisStrata,
            toSize: toSize
        ]
    ]
    withLog(cmd, logMap)

    stub:
    downsampledPairs = "${id}.downsampled.pairs.tsv"
    stub = "touch '${downsampledPairs}'"
    conjunctsArg = conjuncts ? "--conjuncts '${conjuncts.join(' ')}'" : ""
    cisStrataArg = cisStrata ? "--cis-strata '${cisStrata.join(' ')}'" : ""
    toSizeArg = toSize ? "--to-size ${toSize}" : ""
    cmd = "hich downsample ${conjunctsArg} ${cisStrataArg} --orig-stats ${statsFrom} --target-stats ${statsTo} ${toSizeArg} ${fullPairs} ${downsampledPairs}"
    logMap = [
        task: "HichDownsample",
        output: downsampledPairs,
        input: [
            id: id,
            fullPairs: fullPairs,
            statsFrom: statsFrom,
            statsTo: statsTo,
            conjuncts: conjuncts,
            cisStrata: cisStrata,
            toSize: toSize
        ]
    ]
    stubLog(stub, cmd, logMap)
}
/*
    Deduplicate pairs files
*/

process PairtoolsDedup {
    publishDir params.general.publish.parse ? params.general.publish.dedup : "results",
               saveAs: {params.general.publish.dedup ? it : null},
               mode: params.general.publish.mode
    
    container params.general.hichContainer
    label 'pairs'

    input:
    tuple val(id), path(pairs), val(singleCell), val(maxMismatch), val(method), val(pairtoolsDedupParams)

    output:
    tuple val(id), path(deduplicated)

    shell:
    if (!pairtoolsDedupParams) pairtoolsDedupParams = []
    if (singleCell) pairtoolsDedupParams += ["--extra-col-pair cellID cellID --backend cython --chunksize 1000000"]
    if (maxMismatch != null) pairtoolsDedupParams += ["--max-mismatch ${maxMismatch}"]
    if (method != null) pairtoolsDedupParams += ["--method ${method}"]
    deduplicated = "${id}.dedup.pairs.gz"

    cmd = "pairtools dedup --output '${deduplicated}' ${pairtoolsDedupParams.join(' ')}  --nproc-in ${task.cpus} --nproc-out ${task.cpus} '${pairs}'"
    logMap = [task: "PairtoolsDedup", output: deduplicated, input: [id: id, pairs: pairs, singleCell: singleCell, maxMismatch: maxMismatch, method: method, pairtoolsDedupParams: pairtoolsDedupParams]]
    withLog(cmd, logMap)

    stub:
    deduplicated = "${id}.dedup.pairs.gz"
    stub = "touch '${deduplicated}'"
    if (!pairtoolsDedupParams) pairtoolsDedupParams = []
    if (singleCell) pairtoolsDedupParams += ["--extra-col-pair cellID cellID --backend cython --chunksize 1000000"]
    if (maxMismatch != null) pairtoolsDedupParams += ["--max-mismatch ${maxMismatch}"]
    if (method != null) pairtoolsDedupParams += ["--method ${method}"]
    deduplicated = "${id}.dedup.pairs.gz"

    cmd = "pairtools dedup --output '${deduplicated}' ${pairtoolsDedupParams.join(' ')}  --nproc-in ${task.cpus} --nproc-out ${task.cpus} '${pairs}'"
    logMap = [task: "PairtoolsDedup", output: deduplicated, input: [id: id, pairs: pairs, singleCell: singleCell, maxMismatch: maxMismatch, method: method, pairtoolsDedupParams: pairtoolsDedupParams]]
    stubLog(stub, cmd, logMap)
}

/*
    Merge pairs files to one level higher (techreps -> bioreps -> conditions)
*/
process PairtoolsMerge {
    container params.general.hichContainer
    label 'pairs'
    
    input:
    tuple val(id), path(toMerge)

    output:
    tuple val(id), path(merged)

    shell:
    merged = "${id}.merged.pairs.gz"
    toMerge = toMerge.collect { "'${it}'" }
    cmd = "pairtools merge --output '${merged}'  --nproc-in ${task.cpus} --nproc-out ${task.cpus} ${toMerge.join(' ')}"
    logMap = [task: "PairtoolsMerge", input: [id: id, toMerge: toMerge], output: [merged: merged]]
    withLog(cmd, logMap)

    stub:
    merged = "${id}.merged.pairs.gz"
    stub = "touch '${merged}'"
    toMerge = toMerge.collect { "'${it}'" }
    cmd = "pairtools merge --output '${merged}'  --nproc-in ${task.cpus} --nproc-out ${task.cpus} ${toMerge.join(' ')}"
    logMap = [task: "PairtoolsMerge", input: [id: id, toMerge: toMerge], output: [merged: merged]]
    stubLog(stub, cmd, logMap)
    
}


/*
    "Raw samples" are those that have not been assigned to any aggregation profile.

    Aggregation profiles are implemented by copying samples and tagging the copy with
    the aggregation profile and level. This is the function of the CreateAggregateSampleProfiles
    workflow.

    Aggregation profiles are defined in the Nextflow params under params.aggregate.

    Take:
        levelSamples: channel[map[str, Any]] -- the samples filtered for the current level (i.e. techreps, bioreps, or conditions)
        levelParams: map[str, str] -- maps general keywords to level-specific keywords under which relevant values will be stored
    
    Emit:
        result: channel[map[str, Any]] -- the samples after the level-specific params have been applied
        levelParams: map[str, str] -- the same levelParams taken
*/
workflow CreateAggregateSampleProfiles {
    take:
    levelSamples
    levelParams

    main:

    // Store current samples as the result in case there are no aggregation profiles defined
    levelSamples | set{result}

    // Only if "aggregate" is defined in params do we aggregate anything
    if (params.containsKey("aggregate")) {
        // Separate out raw samples
        levelSamples | branch {
            YES: it.aggregateProfileName == null
            NO: true
        } | set{raw}

        /*
            Create a copy of each raw sample for each aggregation profile

            First, reshape the profiles to a columnar format like [profileName: [list of profile names], profileParams: [list of profile param maps]]
            Then convert to a row format and integrate with the samples, ultimately resulting in samples with new parameters:
            [original sample attributes] + [aggregateProfileName: profileName, profileParams: profileParams]

            We give each newly created sample profile-specific id.
        */
        aggregateProfiles = channel.of(params.aggregate)
        // Convert params.aggregate to columnar format with columns "profileName" and "profileParams"
        // profileParams is a map[str, Any] of parameters determining how the sample will be aggregated
        | map{[profileName: it.keySet().toList(), profileParams: it.values().toList()]}

        // Convert the aggregateProfile parameters, which are now in columnar format,
        // to a rows format. Then convert from a single list of maps to a channel of map elements.
        // Then take the Cartesian product of the aggregate profile parameters with the
        // raw samples for the current aggregation level.
        | map{rows(it)}
        | flatten
        | combine(raw.YES)

        // Separate out the profile and sample for each combination and add the profile's parameters
        // to the sample parameters
        | map {
            profile = it[0]
            sample = it[1]
            sample += [aggregateProfileName: profile.profileName] + profile.profileParams
            sample
        }
        | map{
            // We put levelParams.downsampleToMeanDistribution in parentheses to make Groovy treat it as a variable
            // for the purposes of setting a map key rather than as a string.

            downsampleToMeanDistributionGroup = it.subMap((levelParams.downsampleToMeanDistribution)) ?: null

            // Set a new ID for the sample. Also set the downsampleToMeanDistribution submap in as well.
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
        YES: (
            // Select samples that have an aggregateProfileName and that include
            // a value of either downsampleToMeanDistribution or downsampletoSize
            // that requires processing.
            it.aggregateProfileName != null 
            && (
                it.get(levelParams.downsampleToMeanDistribution) != null
                || !(it.get(levelParams.downsampleToSize) in [1.0, null])
            )
        )
        NO: true
    } | set{downsample}

    // Compute the number of reads downsampled from
    downsample.YES
        | map{tuple(it.id, it.latestPairs, it.get(levelParams.readConjuncts), it.get(levelParams.cisStrata))}
        | HichStats
        | map{[id:it[0], (levelParams.downsampleStatsFrom):it[1]]}
        | set{downsampleStatsFromResult}

    // Join the results to the sample map
    pack(downsample.YES, downsampleStatsFromResult) | set {fromStatsCalculated}

    // Compute the number of reads downsampled to
    // These are converted to columnar format as input to HichStatsAggregate,
    // then converted back to rows format and flattened into channel elements
    groupHashMap(fromStatsCalculated, [levelParams.downsampleToMeanDistribution, 'aggregateProfileName'])
        | map{columns(it, ["dropNull":true])}
        | map{tuple(it.id, it.get(levelParams.downsampleStatsFrom), it.outlier)}
        | HichStatsAggregate
        | map{[id:it[0], (levelParams.downsampleStatsTo):it[1]]}
        | map{rows(it)}
        | flatten
        | set {downsampleToGroupMinResult}

    // Join the results to the sample map
    pack(fromStatsCalculated, downsampleToGroupMinResult) | set{toGroupMinResult}

    // Downsample from the number of reads in each bin in levelParams.downsampleStatsFrom to levelParams.downsampleStatsTo
    toGroupMinResult
        | map{
            tuple(it.id, it.latestPairs, it.get(levelParams.downsampleStatsFrom), it.get(levelParams.downsampleStatsTo),
              it.get(levelParams.readConjuncts), it.get(levelParams.cisStrata), it.get(levelParams.downsampleToSize))}
        | HichDownsample
        | map{[id:it[0], (levelParams.downsamplePairs):it[1], latest:it[1], latestPairs:it[1]]}
        | set{downsampledResult}

    // Join the results to the sample map
    pack(toGroupMinResult, downsampledResult) | set {downsampled}

    // Add back in the non-downsampled samples
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

    // Get samples not specifically excluded from merging and where levelParams.doMerge is true
    levelSamples | branch {
        YES: it.includeInMerge != false && it.get(levelParams.doMerge)
        NO: true}
    | set{merge}

    // Merge the pairs files.
    groupHashMap(merge.YES, levelParams.mergeGroupIdentifiers)
        | map{coalesce(columns(it, ["dropNull":true]))}
        | map{tuple(constructIdentifier(coalesce(it, "_drop")), it.latestPairs)}
        | PairtoolsMerge
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{fromMerge}

    // Group the merged result by the mergeGroupIdentifiers, then coalesce common values
    // to a single value, dropping any null or heterogeneous values. The other
    // common values are kept as inherited merge attributes. Then add an ID
    // for the merge.
    groupHashMap(merge.YES, levelParams.mergeGroupIdentifiers)
        | map{coalesce(columns(it, ["dropNull":true]), '_drop')}
        | map{it += [id:constructIdentifier(it)]}
        | set{inheritedMergeAttributes}
    
    // Combine the values from the merge process itself and the inherited merge attributes
    pack(inheritedMergeAttributes, fromMerge) | set{merged}

    // Add to the set of samples
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

    // Extract samples that aren't already duplicated and where they are marked for deduplication
    levelSamples
        | branch {
            YES: !it.alreadyDeduplicated && levelParams.levelFilter(it) && it.get(levelParams.doDedup)
            NO: true
    } | set{deduplicate}

    // Deduplicate the samples and get the resulting attributes
    deduplicate.YES
        | map{tuple(it.id, it.latestPairs, it.dedupSingleCell, it.dedupMaxMismatch, it.dedupMethod, it.pairtoolsDedupParams)}
        | PairtoolsDedup
        | map{[id:it[0], dedupPairs:it[1], latest:it[1], latestPairs:it[1], alreadyDeduplicated:true]}
        | set{dedupResult}

    // Combine with the original sample
    pack(deduplicate.YES, dedupResult) | concat(deduplicate.NO) | set{result}

    emit:
    result
}

/*
Context:

| Select | AggregateTechreps | AggregateBioreps | AggregateConditions | HicMatrix

*/

workflow AggregateTechreps {
    take:
    samples

    main:

    // Create techreps-specific versions of the generic parameters
    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'techrepDedupMethod',
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

    // Separate the techrep-level samples where we are to aggregate the techreps and that have pairs
    samples | branch {
        techrep: !skip("aggregate") && !skip("aggregateTechreps") && isTechrep(it) && (it.pairs || it.latestPairs)
        other: true
    } | set{sampleType}

    // Downsample, merge and deduplicate as necessary
    CreateAggregateSampleProfiles(sampleType.techrep, aggregateParams)
    | DownsamplePairs
    | MergePairs
    | DeduplicatePairs
    | concat(sampleType.other)
    | set{samples}

    // Early stopping
    samples = emptyOnLastStep("aggregateTechreps", samples)

    emit:
    samples
}

workflow AggregateBioreps {
    take:
    samples

    main:

    // Separate the biorep-level samples where we are to aggregate the bioreps and that have pairs
    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'biorepDedupMethod',
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

    // Separate the bioreps-level samples where we are to aggregate the bioreps and that have pairs
    samples | branch {
        biorep: !skip("aggregate") && !skip("aggregateBioreps") && isBiorep(it) && (it.pairs || it.latestPairs)
        other: true
    } | set{sampleType}

    // Downsample, merge and deduplicate as necessary
    CreateAggregateSampleProfiles(sampleType.biorep, aggregateParams)
    | DownsamplePairs
    | MergePairs
    | DeduplicatePairs
    | concat(sampleType.other)
    | set{samples}

    // Early stopping
    samples = emptyOnLastStep("aggregateBioreps", samples)

    emit:
    samples
}

workflow AggregateConditions {
    take:
    samples

    main:

    // Separate the condition-level samples where we are to aggregate the bioreps and that have pairs
    aggregateParams = [
        dedupSingleCell: 'dedupSingleCell',
        dedupMaxMismatch: 'dedupMaxMismatch',
        dedupMethod: 'conditionDedupMethod',
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

    // Separate the condition-level samples where we are to aggregate the bioreps and that have pairs
    samples | branch {
        condition: !skip("aggregate") && !skip("aggregateConditions") && isCondition(it) && (it.pairs || it.latestPairs)
        other: true
    } | set{sampleType}
    
    // Downsample, merge and deduplicate as necessary
    CreateAggregateSampleProfiles(sampleType.condition, aggregateParams)
    | DownsamplePairs
    | MergePairs
    | DeduplicatePairs
    | concat(sampleType.other)
    | set{samples}

    // Early stopping
    samples = emptyOnLastStep("aggregateConditions", samples)

    emit:
    samples
}
