include {emptyOnLastStep; skip} from '../util/cli.nf'
include {keyUpdate} from '../util/keyUpdate.nf'
include {isExistingFile} from '../util/files.nf'
include {withLog; stubLog} from '../util/logs.nf'

process StageGenomeReference {
    /*  When a URL is passed to a Nextflow function, the resource will be
        automatically downloaded and staged by Nextflow. This is a dummy
        function used to download a unique reference and intentionally has
        no content.
    */
    publishDir params.general.publish.genomeReference ? params.general.publish.genomeReference : "results",
               saveAs: {params.general.publish.genomeReference ? it : null},
               mode: params.general.publish.mode
    label 'smallResource'

    input:
    tuple val(assembly), path(uri)

    output:
    tuple val(assembly), path(uri)

    shell:
    cmd = ":"
    logMap = [task: "StageGenomeReference", input: [assembly: assembly, uri: uri]]
    withLog(cmd, logMap)
}

workflow GenomeReference {
    /*
        We will mainly aim to download resources from NCBI:

        https://www.ncbi.nlm.nih.gov/datasets/

        See:
        https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
    */
    take:
        samples

    main:
        /*  
            1. Convert 'assembly' to url key if no existing reference provided
            2. Reference (URL or local file) cast to file(reference) which
               causes Nextflow to automatically stage URLs for download and
               then treat them as local files going forward.
        */
        synonyms = [
         "t2t-hcm13v2.0":"t2t-hcm13v2.0",
         "t2t-chm13":"t2t-hcm13v2.0",
         "t2t":"t2t-hcm13v2.0",
         "hg38":"hg38",
         "GRCh38":"hg38",
         "mm10":"mm10",
         "dm6":"dm6",
         "galGal5":"galGal5",
         "bGalGal5":"galGal5",
         "danRer11":"danRer11",
         "ce10":"ce10"]

        /*
            Reference genome URLs for human and for common model organisms
        */
        urls = ["hg38":"https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz",
                "mm10":"https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz",
                "dm6":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/215/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.gz",
                "galGal5":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/027/408/225/GCA_027408225.1_bGalGal5.pri/GCA_027408225.1_bGalGal5.pri_genomic.fna.gz",
                "danRer11":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/035/GCF_000002035.6_GRCz11/GCF_000002035.6_GRCz11_genomic.fna.gz",
                "ce10":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/985/GCF_000002985.6_WBcel235/GCF_000002985.6_WBcel235_genomic.fna.gz",
                "t2t-hcm13v2.0":"https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/914/755/GCF_009914755.1_T2T-CHM13v2.0/GCF_009914755.1_T2T-CHM13v2.0_genomic.fna.gz"]
        
        /*
            1. Filter for samples tagged with datatype 'fastq'
            2. Group by 'assembly' and 'reference'
            3. Filter for unspecified or non-existent unique files/URLs
            4. Get URL for missing references
            5. Download
            6. Set file as output path
        */
    samples
        | filter{!skip("genomeReference") && !isExistingFile(it.genomeReference)}
        | map{
            errorMessage = """
            Sample ${it.id} with assembly '${it.assembly}' had genomeReference '${it.genomeReference}' which is either null or nonexistent.
            Hich can automatically download genomeReference for common model organisms based on assembly nickname, but this only supports
            the following options: ${synonyms.keySet()}. You can search on https://www.ncbi.nlm.nih.gov/home/genomes/ for genomes for your organism,
            manually download, and specify the path to the filename in the sample file under the genomeReference column. 
            """

            it.genomeReference = it.genomeReference ?: urls.get(synonyms.get(it.assembly))
            
            assert it.genomeReference, errorMessage
            it
        }
        | map{tuple(it.assembly, it.genomeReference)}
        | unique
        | StageGenomeReference
        | map{assembly, genomeReference -> [assembly: assembly, genomeReference: genomeReference]}
        | set{result}

    keyUpdate(samples, result, "assembly") | set{samples}

    samples = emptyOnLastStep("genomeReference", samples)

    emit:
        samples
}