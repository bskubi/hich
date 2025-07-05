include {createCompositeStrategy; filterSamplesByStrategy} from '../util/analysisPlans.nf'
include {skip} from '../util/cli.nf'
include {withLog; stubLog} from '../util/logs.nf'

process CallCompartments {
    publishDir "results/compartments",
               mode: params.general.publish.mode

    label 'features'


    input:
    tuple val(id), path(genomeReference), path(matrix), val(resolution), val(hichCompartmentsParams), val(n_eigs)

    output:
    tuple val(id), path("*.cis.bw"), path("*.cis.vecs.tsv"), path("*.cis.lam.txt"), path("*.phase.bed")

    shell:
    n_eigs = n_eigs ? n_eigs : 3
    cmd = "hich fasta gc ${genomeReference} ${resolution} > ${id}.phase.bed && cooltools eigs-cis --phasing-track ${id}.phase.bed --n-eigs ${n_eigs} --out-prefix ${id}_compartments --bigwig ${matrix}::/resolutions/${resolution}"
    
    logMap = [task: "CallCompartments", input: [id: id, genomeReference: genomeReference, matrix: matrix, resolution: resolution, hichCompartmentsParams: hichCompartmentsParams], output: [cisBW: "${id}.cis.bw", cisVecs: "${id}.cis.vecs.tsv", cisLam: "${id}.cis.lam.txt", phase: "${id}.phase.bed"]]
    withLog(cmd, logMap)

    stub:
    n_eigs = n_eigs ? n_eigs : 3
    cmd = "hich fasta gc ${genomeReference} ${resolution} > ${id}.phase.bed && cooltools eigs-cis --phasing-track ${id}.phase.bed --n-eigs ${n_eigs} --out-prefix ${id}_compartments --bigwig ${matrix}::/resolutions/${resolution}"
    stub = "touch ${id}_compartments.bw ${id}.cis.bw ${id}.cis.vecs.tsv ${id}.cis.lam.txt ${id}.phase.bed"
    logMap = [task: "CallCompartments", input: [id: id, genomeReference: genomeReference, matrix: matrix, resolution: resolution, hichCompartmentsParams: hichCompartmentsParams], output: [cisBW: "${id}.cis.bw", cisVecs: "${id}.cis.vecs.tsv", cisLam: "${id}.cis.lam.txt", phase: "${id}.phase.bed"]]
    withLog(cmd, logMap, stub)
}

workflow CompartmentScores {
    take:
    samples

    main:
    

    if (!skip("compartments")) {
        params.compartments.each {
            planName, analysisPlan ->

            strategy = createCompositeStrategy(analysisPlan.sampleSelectionStrategy, params.sampleSelectionStrategies)

            filterSamplesByStrategy(samples, strategy)
                | map{
                    sample ->
                    tuple(sample.id, sample.genomeReference, sample.latestMatrix, analysisPlan.resolution, analysisPlan.hichCompartmentsParams, analysisPlan.n_eigs)
                }
                | CallCompartments
        }
    }

    emit:
    samples
}
