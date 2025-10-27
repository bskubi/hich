workflow groupSourcesByTarget {
    take:
    ch_sources  // Channel: [source_id, source_record_tuple] - source_id must be unique
    ch_links    // Channel: [source_id, target_id] - source_id can be duplicated
    
    main:
    /*
      Problem: Group source records by target ID based on links.
      `join` requires unique keys in both channels (fails on ch_links).
      Workaround: Use `cross` for one-to-many mapping based on unique source_id,
      then map to re-key by target_id, and finally groupTuple.
     */
    ch_sources
        | cross( ch_links ) // -> [source_id, source_record_tuple, target_id]
        | map { sources, links ->
            // Re-key by target_id, pass through the *entire* original source record tuple
            [links[1], sources]
        }   // -> [target_id, source_record_tuple]
        | groupTuple // -> [target_id, list_of_source_record_tuples]
        | set { ch_grouped_data }

    emit:
    ch_grouped_data // Channel: [target_id, list_of_source_record_tuples]
}