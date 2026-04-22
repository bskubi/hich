include { BWAMEM } from '../modules/local/bwamem/main.nf'
include { BWAMEM2 } from '../modules/local/bwamem2/main.nf'
include { BWAMETH } from '../modules/local/bwameth/main.nf'
include { MERGE_DEDUP } from '../modules/local/mergeDedup/main.nf'
include { EXTRACT_DNA_METHYLATION } from '../modules/local/extractDNAMethylation/main.nf'
include { EXTRACT_3D_PAIRS } from '../modules/local/extract3DPairs/main.nf'

def to_bool(map, key) {
    falsey = [null, "false", "0", "no", "", "n"]
    return [(key): !(map.get(key).toLowerCase() in falsey)]
}

workflow Hich {
    channel.fromPath(params.sampleFile)
        | splitCsv(header: true, sep: "\t") 
        | map{
            sample ->
            // Validate ID
            sample += [ID: sample.ID.trim()]
            if (!sample.ID) { error("Sample ID missing or is only whitespace for ID ${sample.ID}")}

            // Validate flags (convert to boolean)
            ["isInterleaved"].each{ flag -> sample += to_bool(sample, flag) }

            // Validate fastqs
            //  a-b. 'fastqs' is list of files from fastq, fastq1, fastq2 or not in map
            //  c. 'fastq', 'fastq1', 'fastq2' are removed from sample map
            //  d. All given fastq files exist'
            sample += [fastqs: [sample.fastq, sample.fastq1, sample.fastq2].flatten().findResults{it}.collect{file(it)}]
            if (!sample) { sample -= sample.subMap("fastqs") }
            sample -= sample.subMap(["fastq", "fastq1", "fastq2"])
            sample.getOrDefault("fastqs", []).each{ if (!it.exists()) {error("${it} does not exist for ID ${sample.ID}")}}
            
            
            // Validate SAM
            //  a-c. 'sam' is list containing one existing file or not in map
            //  d. 'fastqs' or 'sam' given but not both
            sample += [sam: [sample.sam].flatten().findResults{it}.collect{file(it)}]
            if (!sample.sam) { sample -= sample.subMap("sam") }
            sample.getOrDefault("sam", []).each{ if (!it.exists()) {error("${it} does not exist for ID ${sample.ID}")}}
            if (sample.fastqs && sample.sam) {error("Sample must have either fastqs or sam/bam/cram, not both, for ID ${sample.ID}")}
            if (!sample.fastqs && !sample.sam) {error("Sample must have fastqs or sam for ID ${sample.ID}")}
            
            // Validate index
            //  a-b. indexDir and indexPfx must both be present or absent
            //  c. If fastqs are given, index must be given
            //  d-f. indexDir is existing directory
            //  g. Required index files are present
            if (!sample.indexDir && sample.indexPfx) { error("Missing indexDir for indexPfx ${sample.indexPfx} for ID ${sample.ID}")}
            if (sample.indexDir && !sample.indexPfx) { error("Missing indexPfx for indexDir ${sample.indexDir} for ID ${sample.ID}")}
            if (sample.fastqs && !(sample.indexDir && sample.indexPfx)) { error("To align fastqs, must include indexDir and indexPfx for ID ${sample.ID}")}
            if (sample.indexDir) { sample += [indexDir: file(sample.indexDir)] }
            if (sample.indexDir && !sample.indexDir.exists()) {error("indexDir ${sample.indexDir} does not exist for ID ${sample.ID}")}
            if (sample.indexDir && !sample.indexDir.isDirectory()) {error("indexDir ${sample.indexDir} is not a directory for ID ${sample.ID}")}
            if (params.containsKey("DNAMethylation")) {
                // Combine prefix, common suffixes, and all required final suffixes (bwameth index-mem2)
                ["", ".amb", ".ann", ".pac", ".bwt", ".sa"].each{ sfx ->
                    path = file(sample.indexDir) / ("${sample.indexPfx}.bwameth.c2t${sfx}")
                    if (!path.exists()) {
                        error("Required index file for bwameth ${path} not found using prefix ${sample.indexPfx}.")
                    }
                }

            } else {
                // Combine prefix, common suffixes, and all required final suffixes (bwa-mem2)
            }

            // Validate dedupGroup
            //  a. Must not be all whitespace
            //  b. dedupGroup key is present only if merge is performed
            if (sample.dedupGroup && !sample.dedupGroup.trim()) { error("Sample dedupGroup is all-whitespace for ID ${sample.ID}")}
            if (sample.containsKey("dedupGroup") && !sample.dedupGroup) { sample -= sample.subMap("dedupGroup")}
            
            sample
        }
        | set {samples}

    // Validate IDs and dedupGroups
    //  a. All IDs are unique
    //  b. No IDs match dedupGroups
    samples | map{it.ID} | set {sampleIDs}
    samples | filter{it.dedupGroup} | map{it.dedupGroup} | set {dedupGroups}
    sampleIDs
        | collect
        | map{ 
            collisions = it.countBy{it}.findAll{it.value > 1}
            if (collisions) {error("ID collision. Collision on IDs: ${collisions.keySet()}")            }
        }
    
    sampleIDs
        | join(dedupGroups)
        | map{ match -> error("ID '${match}' matched dedupGroup '${match}'. IDs must not conflict with dedupGroups.")}

    // Use ID as dedupGroup for samples lacking a dedupGroup.
    // Valid since IDs and dedupGroups are guaranteed to not match.
    samples | map { if (!it.dedupGroup) { it += [dedupGroup: it.ID] }; it } | set{samples}

    // Validate how to align and extract attributes
    //  a. Ensure data is explicitly labeled '3DPairs' and/or 'DNAMethylation'
    //  b. Use BWAMETH if --DNAMethylation specified, BWAMEM2 otherwise.
    //  c. Use --3DPairs flag to align chimeric reads
    is3DPairs = params.containsKey("3DPairs")
    isDNAMethylation = params.containsKey("DNAMethylation")
    if (!is3DPairs && !isDNAMethylation) { error("Must run Hich with --3DPairs and/or --DNAMethylation.")}

    ALIGN = isDNAMethylation ? BWAMETH : (params.containsKey("useBwaMem") ? BWAMEM : BWAMEM2)
    samples
        | filter{ it.fastqs }
        | map{ tuple(it.ID, it.fastqs, it.indexDir, it.indexPfx, is3DPairs, it.isInterleaved) }
        | ALIGN
        | set{alignOutput}  // Cannot use ALIGN.output

    // Extract to (ID, sam) tuple matching format of aligner process output
    samples | filter{ it.sam } | map{tuple(it.ID, it.sam)} | set{samInput}

    // Mix in results
    //  a. Combine results to single channel of (ID, sam) tuples
    //  b. Join with samples on ID
    alignOutput | mix(samInput) | set {samData}

    samples
        | map{tuple(it.ID, it)}
        | join(samData)
        | map{ID, sample, sam -> sample + [sam: sam]}
        | set {samples}

    // Merge and deduplicate SAM files
    // Note 1: All samples now have a dedupGroup equal to ID if not provided by user.
    // Note 2: All samples in a dedupGroup have identical values required for downstream analysis
    configEditSAM = file(params.get("configEditSAM"))

    samples
        | map{tuple(it.dedupGroup, it)}
        | groupTuple
        | map {dedupGroup, samples ->
            tuple(dedupGroup, samples.collect{it.sam}.flatten(), samples[0].dedupTag, configEditSAM)
        }
        | MERGE_DEDUP
        | set{mergeDedupOutput} // Can use MERGE_DEDUP.output but use this for consistency


    // Extract DNA methylation
    if (isDNAMethylation && !params.containsKey("skipExtractDNAMethylation")) {
        if (!params.containsKey("configExtractDNAMethylation")) { error("Must include config file with --configExtractDNAMethylation")}
        if (!params.containsKey("genomeReference")) { error("Must include FASTA genome reference with --genomeReference")}
        configExtractDNAMethylation = file(params.get("configExtractDNAMethylation"))
        genomeReference = file(params.get("genomeReference"))
        mergeDedupOutput
            | map{ tuple(*it, genomeReference, configExtractDNAMethylation) }
            | EXTRACT_DNA_METHYLATION
    }

    // // Extract 3D genome
    if (is3DPairs && !params.containsKey("skipExtract3DPairs")) {
        if (!params.containsKey("configExtract3DPairs")) { error("Must include config file with --configExtract3DPairs")}
        if (!params.containsKey("chromsizes")) { error("Must include chromsizes file with --chromsizes")}
        configExtract3DPairs = file(params.get("configExtract3DPairs"))
        chromsizes = file(params.get("chromsizes"))
        addColumns = params.get("addColumns","")
        mergeDedupOutput
            | map { tuple(*it, configExtract3DPairs, chromsizes, addColumns) }
            | EXTRACT_3D_PAIRS
    }
}