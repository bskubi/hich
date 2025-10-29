include { groupSourcesByTarget } from './groupSourcesByTarget.nf'


workflow getSources {
    take:
    ch_targets
    ch_sources

    main:   
    /** Create [source_id, target_id] channel
    */
    ch_targets
        | map { target -> 
            target.source_ids.collect{ source_id ->
                [source_id, target.id]
            }
         }
        | flatMap
        | set { ch_source_target }

    /** Get [target_id, [source_data]] channel from
        ch_source_pairs and ch_source_target.
    */
    groupSourcesByTarget( 
        ch_sources, 
        ch_source_target
    )
        | map{ target_id, sources -> 
            [target_id, sources.sort()]
        }
        | set { ch_merge_output }

    ch_targets
        | map { [it.id, it.source_ids]}
        | join( ch_merge_output )
        | map { target_id, required_source_ids, found_sources -> 
            found_source_ids = found_sources.collect{source_id, data -> source_id}
            found_source_data = found_sources.collect{source_id, data -> data}
            missing = required_source_ids - found_source_ids
            extra = found_source_ids - required_source_ids
            assert !missing && !extra, "Mismatch between required and expected source ids. Required: ${required_source_ids} Found: ${found_source_ids} Missing: ${missing} Extra: ${extra}"
            [target_id, found_source_data.sort()]
        }
        | set { ch_targets }
    

    emit:
    ch_targets
}