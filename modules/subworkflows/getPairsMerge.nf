include { groupSourcesByTarget } from './groupSourcesByTarget.nf'


workflow getPairsMerge {
    take:
    ch_entrypoints
    ch_process_output

    main:
    /** Branch apart the source and target channels from this entrypoint
    */
    ch_entrypoints
        | branch {
            target: it.merge
            source: true
        }
        | set{ ch_entrypoints_merge }
    
    /** Create [source_id, pairs] channel combining entrypoint with
        previous process outputs.
    */
    ch_entrypoints_merge.source
        | map { [it.id, it.pairs] }
        | concat( ch_process_output )
        | set{ ch_source_pairs }
    
    /** Create [source_id, target_id] channel
    */
    ch_entrypoints_merge.target
        | map { target -> 
            target.merge.collect{ source_id ->
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
            [target_id, sources.collect{ id, pairs -> pairs }.sort()]
        }
        | set { ch_pairs_merge }

    emit:
    ch_pairs_merge
}