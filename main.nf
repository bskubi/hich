import groovy.yaml.YamlSlurper
import groovy.json.JsonBuilder
import java.nio.file.Path

include { FASTQ_ALIGN } from './modules/local/fastq/fastq_align/main.nf'
include { BAM_PARSE_PAIRS } from './modules/local/bam/bam_parse_pairs/main.nf'
include { PAIRS_LABEL } from './modules/local/pairs/pairs_label/main.nf'
include { PAIRS_SELECT } from './modules/local/pairs/pairs_select/main.nf'
include { PAIRS_MERGE as PAIRS_MERGE_BEFORE_DEDUP;
          PAIRS_MERGE as PAIRS_MERGE_AFTER_DEDUP  
        } from './modules/local/pairs/pairs_merge/main.nf'
include { PAIRS_DEDUP } from './modules/local/pairs/pairs_dedup/main.nf'
include { PAIRS_COOL_BIN } from './modules/local/pairs/pairs_cool_bin/main.nf'
include { COOL_COARSEN_ADDNORM } from './modules/local/matrix/cool_coarsen_addnorm/main.nf'
include { PAIRS_HIC_BIN_COARSEN_ADDNORM } from './modules/local/pairs/pairs_hic_bin_coarsen_addnorm/main.nf'
include { groupSourcesByTarget } from './modules/subworkflows/groupSourcesByTarget.nf'
include { getPairsMerge } from './modules/subworkflows/getPairsMerge.nf'

/** Respect nf-core naming conventions.
    https://nf-co.re/docs/guidelines/components/subworkflows
*/
workflow {
    manifest_path = params.containsKey("manifest") ? params.manifest : "manifest.yaml"
    manifest_file = file(manifest_path)
    manifest = new YamlSlurper().parse(manifest_file)


    /** Validate hich_version
        Enables us to warn/error if user's running a potentially buggy
        config from earlier versions of Hich.
    */
    accept_versions = ["unstable"]
    if (!manifest.containsKey("hich_version")) {
        error("Manifest ${manifest_file.getAbsolutePath()} has no hich_version field.")
    } else if (!(manifest.hich_version in accept_versions)) {
        error("Unknown version: '${manifest.hich_version}'. Valid choices: ${accept_versions}")
    }


    /** Build all records.
    */
    REMOVE_KEYS = ["hich_version"]
    records = manifest
                .findAll{ k, v -> !(k in REMOVE_KEYS) }
                .collect{ id, config -> [id: id, *: config] }

    channel.fromList(records) | set { ch_all_records }

    /** Split into subchannels based on where the record enters the pipeline.
    */
    VALID_ENTRYPOINTS = [
        "FASTQ_ALIGN", 
        "BAM_PARSE_PAIRS", 
        "PAIRS_LABEL",
        "PAIRS_SELECT",
        "PAIRS_MERGE_BEFORE_DEDUP",
        "PAIRS_DEDUP",
        "PAIRS_MERGE_AFTER_DEDUP",
        "PAIRS_BIN_COARSEN_ADDNORM",
        "COOL_COARSEN_ADDNORM"
    ]
    ch_all_records
        | branch {
            FASTQ_ALIGN: it.entrypoint == "FASTQ_ALIGN"
            BAM_PARSE_PAIRS: it.entrypoint == "BAM_PARSE_PAIRS"
            PAIRS_LABEL: it.entrypoint == "PAIRS_LABEL"
            PAIRS_SELECT: it.entrypoint == "PAIRS_SELECT"
            PAIRS_MERGE_BEFORE_DEDUP: it.entrypoint == "PAIRS_MERGE_BEFORE_DEDUP"
            PAIRS_DEDUP: it.entrypoint == "PAIRS_DEDUP"
            PAIRS_MERGE_AFTER_DEDUP: it.entrypoint == "PAIRS_MERGE_AFTER_DEDUP"
            PAIRS_BIN_COARSEN_ADDNORM: it.entrypoint == "PAIRS_BIN_COARSEN_ADDNORM"
            COOL_COARSEN_ADDNORM: it.entrypoint == "COOL_COARSEN_ADDNORM"
            error: true
        }
        | set { ch_entrypoints }

    ch_entrypoints.error
        | map{ 
            error("No or invalid entrypoint defined for record. Entrypoint: ${it.entrypoint}. Valid entrypoints: ${VALID_ENTRYPOINTS}.\n${it}") 
        }
    
    /** Align FASTQ files
    */
    ch_entrypoints.FASTQ_ALIGN
        | map { [it.id, [file(it.fastq1)], it.fastq2 ? file(it.fastq2) : [], file(it.aligner_index_dir), it.config_fastq_align] }
        | FASTQ_ALIGN

    /** Ingest BAM files
    */
    
    ch_entrypoints.BAM_PARSE_PAIRS
        | map { [it.id, it.bam] }
        | concat( FASTQ_ALIGN.out.bam )
        | set { ch_data_bam_parse_pairs }

    ch_all_records
        | map { [it.id, it.chromsizes, it.config_bam_parse_pairs] }
        | set { ch_config_bam_parse_pairs }
    
    ch_data_bam_parse_pairs
        | join( ch_config_bam_parse_pairs )
        | map{ id, bam, chromsizes, config -> [id, file(bam), file(chromsizes), config]}
        | BAM_PARSE_PAIRS

    /** Label pairs as criteria for selection by PAIRS_SELECT, PAIRS_DEDUP
    */
    
    ch_entrypoints.PAIRS_LABEL
        | map { [it.id, file(it.pairs)]}
        | concat( BAM_PARSE_PAIRS.out.pairs )
        | set { ch_data_pairs_label }
    
    ch_all_records
        | map { [it.id, it.fragment_index, it.config_pairs_label] }
        | set { ch_config_pairs_label }
    
    ch_data_pairs_label
        | join( ch_config_pairs_label )
        | set { ch_pairs_label_input }
    
    ch_pairs_label_input
        | filter { !it[-1].skip }
        | map { id, input_pairs, fragment_index, config_pairs_label ->
            [id, file(input_pairs), file(fragment_index), config_pairs_label]
        }
        | PAIRS_LABEL

    ch_pairs_label_input
        | filter{  it[-1].skip }
        | map { id, pairs, fragment_index, config -> [id, file(pairs)] }
        | concat(PAIRS_LABEL.out.pairs)
        | set { ch_pairs_label }

    // /** Select reads based on traits of individual reads in isolation
    // */

    ch_entrypoints.PAIRS_SELECT
        | map { [it.id, file(it.pairs)] }
        | concat ( ch_pairs_label )
        | set { ch_data_pairs_select }
    
    ch_all_records
        | map { [it.id, it.config_pairs_select] }
        | set { ch_config_pairs_select }
    
    ch_data_pairs_select
        | join ( ch_config_pairs_select )
        | set { ch_pairs_select_input }
    
    ch_pairs_select_input
        | filter { !it[-1].skip }
        | PAIRS_SELECT
    
    ch_pairs_select_input
        | filter { it[-1].skip }
        | map { id, pairs, config -> [id, file(pairs)]}
        | concat(PAIRS_SELECT.out.pairs)
        | set{ ch_pairs_select }

    /** Merge reads before dedup
        Need to collect all reads in the list of read IDs to join.
    */

    getPairsMerge(
        ch_entrypoints.PAIRS_MERGE_BEFORE_DEDUP,
        ch_pairs_select
    )
        | set { ch_data_pairs_merge_before_dedup}

    /** Get merge config channel
    */
    ch_all_records
        | map { [it.id, it.config_pairs_merge] }
        | set { ch_config_pairs_merge_before_dedup }
    
    /** Join in config and merge
    */
    ch_data_pairs_merge_before_dedup
        | join( ch_config_pairs_merge_before_dedup )
        | map { id, pairs, config_pairs_merge -> [id, pairs.collect{file(it)}, config_pairs_merge]}
        | PAIRS_MERGE_BEFORE_DEDUP

    /** Deduplicate
    */

    ch_entrypoints.PAIRS_DEDUP
        | map{ [it.id, it.pairs] }
        | concat(PAIRS_SELECT.out.pairs)
        | concat(PAIRS_MERGE_BEFORE_DEDUP.out.pairs)
        | set { ch_data_dedup }
    
    ch_all_records
        | map { [it.id, it.config_pairs_dedup] }
        | set { ch_config_pairs_dedup }
    
    ch_data_dedup
        | join ( ch_config_pairs_dedup )
        | set { ch_pairs_dedup_input }
    
    ch_pairs_dedup_input
        | filter { !it[-1].skip }
        | map{ id, pairs, config_pairs_dedup -> [id, file(pairs), config_pairs_dedup] }
        | PAIRS_DEDUP
    
    ch_pairs_dedup_input
        | filter { it[-1].skip }
        | map{ id, pairs, config -> [id, file(pairs)]}
        | concat(PAIRS_DEDUP.out.pairs)
        | set{ ch_pairs_dedup }

    /** Merge reads after dedup
    */
    
    getPairsMerge(
        ch_entrypoints.PAIRS_MERGE_AFTER_DEDUP,
        ch_pairs_dedup
    )
        | set { ch_data_pairs_merge_after_dedup}

    /** Get merge config channel
    */
    ch_all_records
        | map { [it.id, it.config_pairs_merge] }
        | set { ch_config_pairs_merge_after_dedup }
    
    /** Join in config and merge
    */
    ch_data_pairs_merge_after_dedup
        | join( ch_config_pairs_merge_after_dedup )
        | map { id, pairs, config_pairs_merge -> [id, pairs.collect{file(it)}, config_pairs_merge]}
        | PAIRS_MERGE_AFTER_DEDUP
    
    /** Create .cool and .mcool contact matrices
    */

    ch_entrypoints.PAIRS_BIN_COARSEN_ADDNORM
        | map { [it.id, it.pairs] }
        | concat( PAIRS_DEDUP.out.pairs )
        | concat( PAIRS_MERGE_AFTER_DEDUP.out.pairs )
        | set { ch_data_pairs_bin_coarsen_addnorm }
    
    ch_all_records
        | map { [it.id, it.chromsizes, it.config_pairs_cool_bin] }
        | set { ch_config_pairs_cool_bin }
    
    ch_data_pairs_bin_coarsen_addnorm
        | join ( ch_config_pairs_cool_bin )
        | filter { !it[-1].skip }
        | map { id, pairs, chromsizes, config_pairs_cool_bin ->
            [id, file(pairs), file(chromsizes), config_pairs_cool_bin]
        }
        | PAIRS_COOL_BIN
    
    ch_entrypoints.COOL_COARSEN_ADDNORM
        | map { [it.id, it.cool] }
        | concat( PAIRS_COOL_BIN.out.cool )
        | set { ch_data_cool_coarsen_addnorm }

    ch_all_records
        | map { [it.id, it.config_cool_coarsen_addnorm] }
        | set { ch_config_coarsen_addnorm }
    
    ch_data_cool_coarsen_addnorm
        | join( ch_config_coarsen_addnorm )
        | filter { !it[-1].skip }
        | COOL_COARSEN_ADDNORM

    /** Create .hic contact matrices
    */

    ch_all_records
        | map { [it.id, it.chromsizes, it.config_pairs_hic_bin_coarsen_addnorm] }
        | set { ch_config_pairs_hic_bin_coarsen_addnorm }
    
    ch_data_pairs_bin_coarsen_addnorm
        | join( ch_config_pairs_hic_bin_coarsen_addnorm )
        | filter { !it[-1].skip }
        | map { id, pairs, chromsizes, config_pairs_hic_bin_coarsen_addnorm ->
            [id, file(pairs), file(chromsizes), config_pairs_hic_bin_coarsen_addnorm]
        }
        | PAIRS_HIC_BIN_COARSEN_ADDNORM
}