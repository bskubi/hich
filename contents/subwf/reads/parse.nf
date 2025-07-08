include {QCReads} from './qcHicReads.nf'
include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {withLog; stubLog} from '../util/logs.nf'

process PairtoolsParse2 {
    publishDir params.general.publish.parse ? params.general.publish.parse : "results",
               saveAs: {params.general.publish.parse ? it : null},
               mode: params.general.publish.mode

    label 'doJobArray'
    label 'pairs'
    conda 'env/cli_env.yml'

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
    
    sortCmd = ["samtools sort -n '${sambam}'"]
    viewCmd = ["samtools view -b '${sambam}'"]

    parse2Cmd = ["pairtools parse2 --flip --assembly '${assembly}' --chroms-path '${chromsizes}' ${parseParams}"]
    sqlCmd = sql ? ["hich pairs sql --memory-limit '${task.memory}' --threads '${task.cpus}' '${sql}' /dev/stdin"] : []
    pairsSortCmd = ["pairtools sort --output '${id}.pairs.gz' --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]

    // Combine the individual commands, then join with a pipe to form the full command
    parseCmd = parse2Cmd + sqlCmd + pairsSortCmd
    sortParse = (sortCmd + parseCmd).join(" | ")
    viewParse = (viewCmd + parseCmd).join(" | ")

    
    cmd = "samtools view '${sambam}' | awk -F'\\t' '{ print \$1 }' | partitioned && ${viewParse} || ${sortParse}"

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
    pairsSortCmd = ["pairtools sort --output '${id}.pairs.gz' --nproc-in ${task.cpus} --nproc-out ${task.cpus}"]

    // Combine the individual commands, then join with a pipe to form the full command
    parseCmd = parse2Cmd + sqlCmd + pairsSortCmd
    sortParse = (sortCmd + parseCmd).join(" | ")
    viewParse = (viewCmd + parseCmd).join(" | ")

    
    cmd = "samtools view '${sambam}' | awk -F'\\t' '{ print \$1 }' | partitioned && ${viewParse} || ${sortParse}"


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

workflow Parse {
    take:
    samples

    main:
    samples
        | filter{!skip("parse") && it.datatype in ["fastq", "sambam"]}
        | map{tuple(it.id, it.sambam, it.chromsizes, it.assembly, it.pairtoolsParse2Params, it.parseSQL, it.subMap("minMapq"))}
        | PairtoolsParse2
        | map{[id:it[0], pairs:it[1], latest:it[1], latestPairs:it[1]]}
        | set{result}
    keyUpdate(samples, result, "id") | set{samples}

    // It might be good to simplify these workflow control steps since they
    // are repeated frequently.
    if ("parse" in params.general.get("qcAfter")) {
        QCReads(samples, "parse")
    }

    samples = emptyOnLastStep("parse", samples)

    emit:
        samples
}
