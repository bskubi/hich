include { groupSourcesByTarget } from './groupSourcesByTarget.nf'


workflow getSources {
    take:
    ch_entrypoints
    ch_process_output
    key_data

    main:
    /** Branch apart the source and target channels from this entrypoint
    */
    ch_entrypoints
        | branch {
            target: it.source_ids
            source: true
        }
        | set{ ch_entrypoints_merge }
    
    /** Create [source_id, pairs] channel combining entrypoint with
        previous process outputs.
    */
    ch_entrypoints_merge.source
        | map { [it.id, it[key_data]] }
        | concat( ch_process_output )
        | set{ ch_source_pairs }
    
    /** Create [source_id, target_id] channel
    */
    ch_entrypoints_merge.target
        | map { target -> 
            target.sources.collect{ source_id ->
                [source_id, target.id]
            }
         }
        | flatMap
        | set { ch_source_target }

    /** Get [target_id, [source_data]] channel from
        ch_source_pairs and ch_source_target.
    */
    groupSourcesByTarget( 
        ch_source_pairs, 
        ch_source_target
    )
        | map{ target_id, sources -> 
            [target_id, sources.sort()]
        }
        | set { ch_pairs_merge_output }

    ch_entrypoints_merge.target
        | map { [it.id, it.sources]}
        | join( ch_pairs_merge_output )
        | map { target_id, required_source_ids, found_sources -> 
            found_source_ids = found_sources.collect{source_id, data -> source_id}
            found_source_data = found_sources.collect{source_id, data -> data}
            missing = required_source_ids - found_source_ids
            extra = found_source_ids - required_source_ids
            assert !missing && !extra, "Mismatch between required and expected source ids. Required: ${required_source_ids} Found: ${found_source_ids} Missing: ${missing} Extra: ${extra}"
            [target_id, found_source_data.sort()]
        }
        | set { ch_pairs_merge }
    

    emit:
    ch_pairs_merge
}