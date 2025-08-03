include {withLog; stubLog} from '../../util/logs.nf'
include {cmdList} from '../../util/dataStructures.nf'

process SPLIT_PAIRS {
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer
    
    input:
    tuple val(id), path(pairs), val(splitColumns), val(sql)

    output:
    tuple val(id), path("__results__/*")

    shell:
    inPattern = "\"__temp__/chrom1={chrom1}/data_0.parquet\""
    outPattern = "\"__results__/CELL_{chrom1}_${id}.cell={chrom1}.pairs.gz\""
    pairsFile = "\"${pairs}\""
    splitColumns = cmdList(splitColumns)
    sql = sql ? "--sql \"${sql}\"" : ""
    
    cmd = "hich pairs partition --in-pattern ${inPattern} --out-pattern ${outPattern} ${sql} ${pairsFile} __temp__ ${splitColumns} && rm -rf __temp__"
    logMap = [task: "SPLIT_PAIRS", output: "__results__/*.pairs.gz", input: [id: id, pairs: pairs, splitColumns: splitColumns, sql: sql]]
    withLog(cmd, logMap)
}