workflow keydiff {
    take:
        left
        right
        by
        how

    main:
    /*
        Motivation and implementation logic:

        When Nextflow channel items are LinkedHashMaps ("hashmaps"), they are
        similar to un-normalized tables. Items are rows, keys represent the
        column name for that row, and values represent cell entries.

        We can therefore do an SQL-like join on them. A left join, for example,
        keeps all 'left' rows, adding 'right' rows where there is a matching]
        key and null or default values where there is no matching 'right' key.

        This workflow is used to determine which key column values are missing
        in the left table, right table, or in neither table.

        'by': the set of hashmap keys (i.e. table columns) to join on
        'keyset': the values in the by keys of a particular channel item (i.e. table row)
        
    */

    /* Extract unique keysets from left and right tables.

    This also adds a 'true' value which will be used in the next step to detect
    unmatched keysets.
    */
    left_keys = left
        | map{it.subMap(by)}
        | unique
        | map{[it, true]}

    right_keys = right
        | map{it.subMap(by)}
        | unique
        | map{[it, true]}
    
    /* Determine left ∩ right, left - right, and right - left keysets
    
    Nextflow join operator matches on first element of each item (the keyset)
    The remainder argument keeps unmatched items, regardless of which
    channel they were in. It adds a null value as a placeholder which is groovy
    falsey. If the keyset is missing in neither (found in both), the output
    will be [keyset, true, true] from the true values added in the previous
    step. If the keyset is missing in the left, [keyset, null, true] will be
    emitted, or [keyset, true, null] if it is missing in the right.
    */
    missing =
        left_keys
        | join(right_keys, remainder: true)
        | branch {
            neither: it[1] &&  it[2] //key present in both
            right:   it[1] && !it[2] //key in left only
            left:   !it[1] &&  it[2] //key in right only
        }

    /*
        Extract the keyset for rows found only in the left or only in the right
        depending on the 'how' argument, so that default values can be added
        where necessary during the join.
    */
    if (how == "right") {
        missing.right | map{it[0]} | set{result}
    } else if (how == "left") {
        missing.left | map{it[0]} | set{result}
    }

    emit:
        result


}

workflow sqljoin {
    take:
    left
    right
    condition

    main:
    /* Similar to a SQL inner join

    Tables are represented as channels emitting LinkedHashMap 'rows'

    Arguments:
    left -- left channel
    right -- right channel
    by -- Key or list of keys to join on
    suffix -- suffix to add to right if necessary to avoid non-key clashes

    Notes:
    1. If there is a conflict (i.e. the left table has a non-key 'data'
    and the right has a non-key 'data' and 'data_right' with the suffix
    '_right', then the right keys with higher index have more suffix copies
    appended to them.

    2. sqljoin takes key type into account, so even if key1 == key2,
    they will not be joined (and will not be treated as conflicting and
    thus distinguished with a suffix) unless their class is the same. One
    way this can be an issue is if their types are GString and String.
    */

    /*
        Motivation and implementation logic:

        When Nextflow channel items are LinkedHashMaps ("hashmaps"), they are
        similar to un-normalized tables. Items are rows, keys represent the
        column name for that row, and values represent cell entries.

        We can therefore do an SQL-like join on them. A left join, for example,
        keeps all 'left' rows, adding 'right' rows where there is a matching]
        key and null or default values where there is no matching 'right' key.
        If two or more 'right' rows match one 'left' row, then multiple
        rows are produced containing one copy of the 'left' items and one
        copy of one of the matching 'right' rows' items.

        Nextflow's 'join' is similar but importantly different from an SQL
        inner join. Other join types are not available. None operate on named
        keys from hashmaps. We implement a closer approximation to the SQL
        full, left, right and inner joins on Nextflow channels containing
        a single LinkedHashMap item here.

        'keyset': the values in the by keys of a particular channel item (i.e. table row)

        The implementation idea is as follows:

        1. Depending on join type, determine missing keysets that would
        inappropriately lead to rows being dropped and add in empty placeholders.
        For example, for a left join all left rows should be kept. So any keysets
        in the left table not present in the right table must be added as
        placeholders to the right table.

        2. Group the left channel by keyset, giving a [keyset: [list of items]]
        channel. This is necessary because the Nextflow join operator returns
        only the first item from the left channel for each specific keyset.

        3. Use Nextflow cross operator to get pairwise combinations of left
        channel items and right channel items with matching keysets. For example,
        if there are two left items and two right items with a particular
        keyset we obtain:
            [[keyset, [left_1, left_2]], [keyset, right_1]]
            [[keyset, [left_1, left_2]], [keyset, right_2]]

        4. Reformat to produce hashmap-formatted elements:
            [left_1 + right_1]
            [left_1 + right_1]
            [left_1 + right_1]
            [left_1 + right_1]

            Note: the matching keyset will be represented once in the
            resulting hashmaps

            Where there are non-keyset matching columsn in the left and right
            tables, we add a suffix to avoid a clash or (if a blank suffix
            is provided) we overwrite.
    */

    by = condition.get("by", null)
    suffix = condition.get("suffix", "_right")
    how = condition.get("how", "left")

    /* 1. Depending on join type, determine missing keysets that would
       inappropriately lead to rows being dropped and add in empty placeholders.
       For example, for a left join all left rows should be kept. So any keysets
       in the left table not present in the right table must be added as
       placeholders to the right table.
    */
    if (how == "full" || how == "left") {
        // Add any keys in left that are missing in right to right.
        missing = keydiff(left, right, by, "right")
        right = right | concat(missing)
    }
    if (how == "full" || how == "right") {
        // Add any keys in left that are missing in right to right.
        missing = keydiff(left, right, by, "left")
        left = left | concat(missing)
    }

    /* 2. Group the left channel by keyset, giving a [keyset: [list of items]]
       channel. This is necessary because the Nextflow join operator returns
       only the first item from the left channel for each specific keyset.
    */

    left = left | map{[it.subMap(by), it]} | groupTuple
    right = right | map{[it.subMap(by), it]}

    /* 3. Use Nextflow cross operator to get pairwise combinations of left
       channel items and right channel items with matching keysets. For example,
       if there are two left items and two right items with a particular
       keyset we obtain:
           [[keyset, [left_1, left_2]], [keyset, right_1]]
           [[keyset, [left_1, left_2]], [keyset, right_2]]
    */

    left
    | cross(right)
    | map{
        /* Reformat item format
           
           from [[keyset, [left_1, left_2]], [keyset, right]]

           to   [left_1 + right]
                [left_2 + right]

            These are hashmaps with keyset represented only once
        */
        result = []

        /* Drop the key (which is still preserved in the hashmaps)
           and extract the hashmaps from the [keyset, hashmap] elements to get:
           left = [leftrow1, leftrow2, ...]
           rightrow = rightrow

           it[0] = [keyset, [left_1, left_2, ...]]
           it[1] = [keyset, right]
           
           drop(1) means to drop 1 item from the beginning of the list, so it
           will drop the keyset.
        */

        /* This is an old implementation. I'm pretty sure it is unnecessary
        to drop and then extract the first element, but I'm leaving it in
        case it introduces bugs.
        */
        //left = it[0].drop(1)[0]       
        //rightrow = it[1].drop(1)[0]
        
        /* The new version just extracts the desired elements directly*/
        left = it[0][1]
        rightrow = it[1][1]

        // Iterate through each leftrow and combine with rightrow
        left.each() {
            leftrow ->
            combined = [:]

            // Put all items (including key) from leftrow into combined
            combined.putAll(leftrow)
            
            // Iterate through each key in rightrow,
            // appending suffix if needed to avoid clashes with leftrow
            // and adding it to the combined row.
            rightrow.each() {
                key, value ->

                // Add non-keyset keys from right to left.
                if (!by.contains(key)) {

                    // Add as many suffixes as necessary to avoid clash
                    while (suffix != "" && combined.containsKey(key)) {
                        key = key + suffix
                    }

                    /* Add the right key value if either criterion holds:
                        - There is no name clash, or
                        - The left table is not declared dominant
                    */
                    if (condition.get("dominant") != "left" || !combined.containsKey(key)) {
                        combined[key] = value
                    }
                }
            }
            // Add the new hashmap to the list of hashmaps
            result.add(combined)
        }

        // List of per-keyset joins
        result
    }
    | collect       // Currently, channel items are lists of per-keyset joins.
    | flatMap{it}   // Collects to a single list of lists, then flattens
    | set{joined}   // Yielding the desired result of single-hashmap join items.
    
    emit:
        joined
}