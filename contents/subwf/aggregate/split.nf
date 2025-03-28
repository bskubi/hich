process splitTechreps {
    input:
    tuple val(id), path(pairs), val(split)

    output:
    tuple val(id), path(downsampledPairs)

    shell:
    "hich pairs partition --in-pattern '__temp__/${split}={${split}}/data_0.parquet' --out-pattern '${id}_${split}.pairs' '${pairs}' __temp__ ${split}"
}

workflow SplitTechreps {
    take:
    samples

    main:

    

    emit:
    samples
}