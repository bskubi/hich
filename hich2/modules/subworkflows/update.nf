workflow update {
    take:
    samples
    new_data

    main:

    samples
        | map{[it.id, it]}
        | join(new_data)
        | map{it[1] + it[2]}
        | set{samples}

    emit:
    samples
}