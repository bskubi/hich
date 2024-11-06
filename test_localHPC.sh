set -euo pipefail

# Different assemblies, enzymes, aligners, aggregation profiles, sample selection strategies, feature calling parameterizations
# Multiple versions for different -profile settings, containerization solutions

# Define the initial part of the command in a variable
base_command="nextflow run contents/main.nf -profile localHPC -c params/apptainer.config -w testwork -resume"

# Tests: Making bwa aligner index and other input resources
$base_command --sampleFile vignettes/akgol2021/one_rep_bwa_mem2.tsv -params-file params/onerep_bulk.yml

# Tests: Tiered merge and deduplication
$base_command --sampleFile vignettes/akgol2021/samples_bwa_mem2.tsv -params-file params/tiered_bulk.yml

# Tests: Using --fastqInterleaved to process interleaved fastq data, setting defaults via command line, BSBolt alignment, BSBolt aligner indexing
$base_command -param-file params/onesample_bulk.yml --fastqInterleaved "vignettes/snmC-seq3/SRR13751346.fastq.gz" --defaults.condition SRR13751346 --defaults.minMapq 10 --defaults.aligner bsbolt --defaults.assembly hg38 --defaults.restrictionEnzymes "NlaIII,MboI" --defaults.bwaFlags "'-SP5Mp'" --lastStep hicMatrix



