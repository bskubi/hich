include {buildCLIOpts} from '../../util/cli.nf'

def buildCmd(planName, mcools, hicrep_combinations_opts, cpus) {
    def output = "${planName}.tsv"
    def logMap = [
        task: "HICREP_COMBINATIONS", 
        input: [
            planName: planName, 
            mcools: mcools, 
            hicrep_combinations_opts: hicrep_combinations_opts
        ],
        output: [hicrep: output]
    ]

    def default_hich_matrix_hicrep_opts = [
        "--n_proc": cpus
    ]
    def hich_matrix_hicrep_opts = hicrep_combinations_opts?.hich_matrix_hicrep_opts ?: [:]
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
    def final_hich_matrix_hicrep_opts = buildCLIOpts(default_hich_matrix_hicrep_opts, hich_matrix_hicrep_opts, remap, null)
    def args = [output] + mcools
    args = args.collect{"'${it}'"}.join(" ")
    def cmd = "hich matrix hicrep ${final_hich_matrix_hicrep_opts} ${args}"

    return [cmd, logMap, output]
}