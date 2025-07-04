include {columns; sortMapList} from './reshape.nf'

workflow GroupToColumnar {
    take:
    chan
    groupBy
    sortBy
    columnsOptions

    main:
    /*
        Group map channel by keys and convert each group to column format

        Example:
            chan channel.of([a: 1, b: 1, id: 2], [a: 1, b: 2, id: 1], [a: 2, b: 3, id: 3])
            groupBy ["a"]
            sortBy ["id"]
            columnsOptions [:]

            returns
                [
                    [a: [1, 1], b: [2, 1], id: [1, 2]],
                    [a: [2], b: [3], id: [3]]
                ]
        
        List groupBy: Keys to group maps by
        List sortBy: Keys to sort within groups for deterministic group outputs
        Map columnsOptions: Map of options to 'columns' method 
    */
    chan
        | map{tuple(it.subMap(groupBy), it)}
        | groupTuple                    
        | map{it[1]}
        | map{ mapList -> sortBy ? sortMapList(mapList, sortBy) : mapList }
        | map{ mapList -> columns(mapList, columnsOptions)}
        | set{chan}

    emit:
    chan
}