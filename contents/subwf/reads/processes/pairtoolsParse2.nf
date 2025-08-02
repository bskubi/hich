include {withLog; stubLog} from '../util/logs.nf'

process PairtoolsParse2 {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null},
               mode: params.general.publish.mode

    label 'pairs'
    tag "$id"
    conda "$projectDir/env/dev_env.yml"
    container params.general.hichContainer

    input:
    tuple val(id), path(sambam), path(chromsizes), val(assembly), val(parseParams), val(sql), val(flags)

    output:
    tuple val(id), path("${id}.pairs.gz")

    shell:
    parseParams = parseParams ?: []

    // Use minMapq as default, but override with manually specified --min-mapq
    if (flags.minMapq instanceof Integer && !parseParams.any{it.contains("--min-mapq")}) {
        parseParams += ["--min-mapq ${flags.minMapq}"]
    } 

    // The parseParams are typically a list of individual pairtools parse flags.
    // Join them separated by spaces to use in the parse2Cmd
    parseParams = parseParams.join(" ")

    // Set up the individual commands in lists to make them easier to combine with pipes into a complete command
    // sambamba is both slower than samtools as of 2017, and also can't pipe to stdout, so we use samtools
    
    sortCmd = ["samtools sort -@ ${task.cpus} -n '${sambam}'"]
    viewCmd = ["samtools view -b '${sambam}'"]

    // Use 2G less than the total memory allocated for the job
    // to a minimum of 2G
    memory = task.memory ? Math.max(task.memory.toGiga() - 2, 2) : 2

    parse2Cmd = ["pairtools parse2 --flip --assembly '${assembly}' --chroms-path '${chromsizes}' ${parseParams}"]
    sqlCmd = sql ? ["hich pairs sql --memory-limit '${memory}' --threads '${task.cpus}' '${sql}' /dev/stdin"] : []
    
    
    pairsSortCmd = ["pairtools sort --output '${id}.pairs.gz' --memory ${memory}G --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]

    // Combine the individual commands, then join with a pipe to form the full command
    parseCmd = parse2Cmd + sqlCmd + pairsSortCmd
    sortParse = (sortCmd + parseCmd).join(" | ")
    viewParse = (viewCmd + parseCmd).join(" | ")
    
    cmd = viewParse

    logMap = [
        task: "PairtoolsParse2",
        input: [
            id: id,
            sambam: sambam,
            chromsizes: chromsizes,
            assembly: assembly,
            parseParams: parseParams,
            sql: sql,
            flags: flags
        ],
        output: [
            pairs: "${id}.pairs.gz"
        ]
    ]

    withLog(cmd, logMap)



    stub:
    stub = "touch '${id}.pairs.gz'"
    parseParams = parseParams ?: []
    sql = sql ?: []

    // Use minMapq as default, but override with manually specified --min-mapq
    if (flags.minMapq instanceof Integer && !parseParams.any{it.contains("--min-mapq")}) {
        parseParams += ["--min-mapq ${flags.minMapq}"]
    } 

    // The parseParams are typically a list of individual pairtools parse flags.
    // Join them separated by spaces to use in the parse2Cmd
    parseParams = parseParams.join(" ")
    sql = sql.join(" ")

    // Set up the individual commands in lists to make them easier to combine with pipes into a complete command
    // sambamba is both slower than samtools as of 2017, and also can't pipe to stdout, so we use samtools
    
    sortCmd = ["samtools sort -n '${sambam}'"]
    viewCmd = ["samtools view -b '${sambam}'"]

    parse2Cmd = ["pairtools parse2 --flip --assembly '${assembly}' --chroms-path '${chromsizes}' ${parseParams}"]
    sqlCmd = sql ? ["hich pairs sql --memory-limit ${} '${sql}' /dev/stdin"] : []
    memory = task.memory ? task.memory.toGiga() : "2"
    pairsSortCmd = ["pairtools sort --output '${id}.pairs.gz' --memory ${memory}G --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]

    // Combine the individual commands, then join with a pipe to form the full command
    parseCmd = parse2Cmd + sqlCmd + pairsSortCmd
    sortParse = (sortCmd + parseCmd).join(" | ")
    viewParse = (viewCmd + parseCmd).join(" | ")

    
    // cmd = "samtools view '${sambam}' | head -n 10000 | awk -F'\\t' '{ print \$1 }' | partitioned && ${viewParse} || ${sortParse}"
    cmd = viewParse


    logMap = [
        task: "PairtoolsParse2",
        input: [
            id: id,
            sambam: sambam,
            chromsizes: chromsizes,
            assembly: assembly,
            parseParams: parseParams,
            sql: sql,
            flags: flags
        ],
        output: [
            pairs: "${id}.pairs.gz"
        ]
    ]

    stubLog(stub, cmd, logMap)
}