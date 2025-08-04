include {formatArg} from '../../util/cli.nf'

def buildCmd(planName, mcools, resolutions, chroms, exclude, chromFilter, h, dBPMax, bDownSample) {
    def output = "${planName}.tsv"
    def resolution = "-r " + resolutions.join(" -r ")

    def cmd = [
        "hich matrix hicrep",
        resolution,
        formatArg("--chrom %s", chroms, ','),
        formatArg("--exclude %s", exclude, ','),
        formatArg("--chrom-filter '%s'", chromFilter, ''),
        formatArg("--h %s", h, ','),
        formatArg("--d-bp-max %s", dBPMax, ','),
        formatArg("--b-downsample %s", bDownSample, ','),
        "'${output}'",
        formatArg("%s", mcools, ' ')
    ]
    cmd = cmd.findAll{it}
    cmd = cmd.join(" ")
    def logMap = [
        task: "HICREP_COMBINATIONS", 
        input: [
            planName: planName, 
            mcools: mcools, 
            resolutions: resolutions, 
            chroms: chroms, 
            exclude: exclude, 
            chromFilter: chromFilter, 
            dBPMax: dBPMax, 
            bDownSample: bDownSample
        ],
        output: [hicrep: output]
    ]
    return [cmd, logMap, output]
}