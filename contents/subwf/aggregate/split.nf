process Split {
    input:
    tuple val(id), path(pairs)

    output:
    tuple val(id), path("*")

    shell:
    inPattern = "\"__temp__/chrom1={chrom1}/data_0.parquet\""
    outPattern = "\"CELL_{chrom1}_${id}.cell={chrom1}.pairs.gz\""
    pairsFile = "\"${pairs}\""

    "hich pairs partition --in-pattern ${inPattern} --out-pattern ${outPattern} ${pairsFile} __temp__ chrom1 && rm -rf __temp__"
}