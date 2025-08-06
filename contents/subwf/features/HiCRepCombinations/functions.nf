include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(planName, mcools, analysisPlan, cpus) {
    def output = "${planName}.tsv"
    def logMap = [
        task: "HICREP_COMBINATIONS", 
        input: [
            planName: planName, 
            mcools: mcools, 
            analysisPlan: analysisPlan
        ],
        output: [hicrep: output]
    ]

    def default_hich_matrix_hicrep_opts = [
        "--n_proc": cpus
    ]
    def hich_matrix_hicrep_opts = analysisPlan?.hich_matrix_hicrep ?: [:]
    def remap = [
        "--resolution": "-r",
        "--h": "-h",
        "--dBPMax": "-m",
        "--d-bp-max": "-m",
        "--bDownSample": "-d",
        "--b-downsample": "-d",
        "--normalization": "-n",
        "--chrom": "-c",
        "--skip-chrom": "-s",
        "--region": "-g",
        "--bed": "-b",
        "--partition": "-p"
    ]
    hich_matrix_hicrep_opts = buildCLIOpts(default_hich_matrix_hicrep_opts, hich_matrix_hicrep_opts, remap, null)
    mcools = mcools.collect{"'${it}'"}.join(" ")
    def cmd = "hich matrix hicrep ${hich_matrix_hicrep_opts} '${output}' ${mcools}"

    return [cmd, logMap, output]
}