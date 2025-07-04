workflow missingKeys {
    take:
    fullSet
    subSet
    submapBy

    main:
    assert submapBy instanceof List, "In submapdiff, submapBy must be of type List, but was type ${submapBy.getClass()} with value ${submapBy}"
    
    /*
        Identify submaps unique to left or right channel

        channel<Map> left
        channel<Map> right
        List by
        String how  
    */

    // Use submapBy to get distinct submaps from left and right channel items
    fullSetKeys = fullSet
        | map{it.subMap(submapBy)}
        | unique
        | map{[it, true]}

    subSetKeys = subSet
        | map{it.subMap(submapBy)}
        | unique
        | map{[it, true]}

    // Identify submaps present in left only or right only
    fullSetKeys
        | join(subSetKeys, remainder: true)
        | filter{it[1] && !it[2]}
        | map{it[0]}
        | set{missing}


    emit:
    missing
}