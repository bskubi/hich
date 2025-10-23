import groovy.yaml.YamlSlurper
import groovy.json.JsonBuilder
import hich.Hich
import java.nio.file.Path

include { FASTQ_ALIGN } from './modules/local/fastq/fastq_align/main.nf'
include { BAM_PARSE_CONTACTS } from './modules/local/bam/bam_parse_contacts/main.nf'
include { PAIRS_PREP } from './modules/local/pairs/pairs_prep/main.nf'
include { PAIRS_LABEL } from './modules/local/pairs/pairs_label/main.nf'
include { PAIRS_SELECT } from './modules/local/pairs/pairs_select/main.nf'
include { PAIRS_MERGE as PAIRS_MERGE_BEFORE_DEDUP;
          PAIRS_MERGE as PAIRS_MERGE_AFTER_DEDUP  
        } from './modules/local/pairs/pairs_merge/main.nf'
include { PAIRS_DEDUP } from './modules/local/pairs/pairs_dedup/main.nf'
include { PAIRS_HIC_BIN_COARSEN_ADDNORM } from './modules/local/pairs/pairs_hic_bin_coarsen_addnorm/main.nf'
include { PAIRS_COOL_BIN } from './modules/local/pairs/pairs_cool_bin/main.nf'
include { COOL_COARSEN_ADDNORM } from './modules/local/matrix/cool_coarsen_addnorm/main.nf'


List<Path> toPathList(List<String> path_str) {
    return path_str.findAll().collect{file(it)}
}

List<Path> toPathList(String path_str) {
    return [file(it)]
}

/** Respect nf-core naming conventions.
    https://nf-co.re/docs/guidelines/components/subworkflows
*/
workflow {
    hich_config_path = params.containsKey("hichConfig") ? params.configFile : "hich_config.yaml"
    hich_config_file = file(hich_config_path)
    hich_config = new YamlSlurper().parse(hich_config_file)
    
    channel.fromList(
        hich_config
        .input_records
        .collect { id, record -> record + [id: id] + hich_config.standard_alignment + hich_config.M129 }
    )
        | set {ch_record_input}
    
    channel.fromList(
        hich_config
        .merge_before_dedup_records
        .collect { id, record -> record + [id: id] + hich_config.standard_alignment + hich_config.M129 }
    )
        | set {ch_record_merge_before_dedup}

    channel.fromList(
        hich_config
        .merge_after_dedup_records
        .collect { id, record -> record + [id: id] + hich_config.standard_alignment + hich_config.M129 }
    )
        | set {ch_record_merge_after_dedup }

    
    ch_record_input
        | branch {
            fastq_input:  it.fastq && !it.sambam && !it.pairs && !it.cool && !it.mcool && !it.hic
            bam_input: (it.sam || it.bam || it.cram) && !it.pairs && !it.cool && !it.mcool && !it.hic
            pairs_input:  it.pairs && !it.cool && !it.mcool && !it.hic
            matrix_input: it.cool || it.mcool || it.hic
            error: true
        }
        | set{ ch_record_input_source }
    
    ch_record_input_source.fastq_input
        | map { [it.id, toPathList(it.fastq)] }
        | set { ch_fastq_input }

    ch_record_input
        | map { [it.id, it.aligner, file(it.aligner_index_dir), it.aligner_index_prefix, it.config_fastq_align] }
        | set { ch_config_fastq_align }
    
    ch_fastq_input
        | join( ch_config_fastq_align )
        | FASTQ_ALIGN

    ch_record_input_source.bam_input
        | map { [it.id, toPathList([it.sam, it.bam, it.cram])] }
        | set { ch_bam_input }

    ch_record_input
        | map { [it.id, it.config_bam_parse_contacts] }
        | set { ch_config_bam_parse_contacts }

    FASTQ_ALIGN.out.bam
        | concat( ch_bam_input )
        | join( ch_config_bam_parse_contacts )
        | BAM_PARSE_CONTACTS

    ch_record_input_source.pairs_input
        | map { [it.id, toPathList(it.pairs)] }
        | set { ch_pairs_input }
    
    ch_record_input
        | map { [it.id, it.config_pairs_prep] }
        | set { ch_config_pairs_prep }

    ch_pairs_input
        | join( ch_config_pairs_prep )
        | PAIRS_PREP

    ch_record_input
        | map { [it.id, it.config_pairs_label] }
        | set { ch_config_pairs_label }
    
    BAM_PARSE_CONTACTS.out.pairs
        | concat(PAIRS_PREP.out.pairs)
        | join( ch_config_pairs_label )
        | PAIRS_LABEL

    ch_record_input
        | map { [it.id, it.config_select] }
        | set { ch_config_pairs_select }
    
    PAIRS_LABEL.out.pairs
        | join( ch_config_pairs_select)
        | PAIRS_SELECT

    ch_record_input
        | map { [it.id, it.config_pairs_merge_before_dedup] }
        | set { ch_config_pairs_merge_before_dedup }

    PAIRS_SELECT.out.pairs
        | join ( ch_config_pairs_merge_before_dedup )
        | map{ [*it, "_before_dedup"] }
        | PAIRS_MERGE_BEFORE_DEDUP

    ch_record_input
        | map { [it.id, it.config_pairs_dedup] }
        | set { ch_config_pairs_dedup }
    
    PAIRS_SELECT.out.pairs
        | concat(PAIRS_MERGE_BEFORE_DEDUP.out.pairs)
        | join ( ch_config_pairs_dedup )
        | PAIRS_DEDUP
    
    ch_record_input
        | map { [it.id, it.config_pairs_merge_after_dedup] }
        | set { ch_config_pairs_merge_after_dedup }

    PAIRS_DEDUP.out.pairs
        | join ( ch_config_pairs_merge_after_dedup )
        | map{ [*it, "_after_dedup"] }
        | PAIRS_MERGE_AFTER_DEDUP
    
    ch_record_input
        | map { [it.id, it.config_pairs_cool_bin] }
        | set { ch_config_pairs_cool_bin }

    PAIRS_DEDUP.out.pairs
        | concat( PAIRS_MERGE_AFTER_DEDUP.out.pairs )
        | set { ch_output_pairs }
    
    ch_output_pairs
        | join ( ch_config_pairs_cool_bin )
        | PAIRS_COOL_BIN

    ch_record_input_source.matrix_input
        | filter{ it.cool && !it.mcool }
        | map { [it.id, toPathList(it.cool)] }
        | set { ch_cool_input }

    ch_record_input
        | map { [it.id, it.config_cool_coarsen_addnorm] }
        | set { ch_config_cool_coarsen_addnorm }
    
    PAIRS_COOL_BIN.out.cool
        | concat( ch_cool_input )
        | join( ch_config_cool_coarsen_addnorm ) 
        | COOL_COARSEN_ADDNORM
    
    ch_record_input
        | map { [it.id, it.config_hic_bin_coarsen_addnorm] }
        | set { ch_config_hic_bin_coarsen_addnorm }
    
    ch_output_pairs
        | join ( ch_config_hic_bin_coarsen_addnorm )
        | PAIRS_HIC_BIN_COARSEN_ADDNORM
}