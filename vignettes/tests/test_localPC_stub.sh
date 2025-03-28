set -euo pipefail

# Different assemblies, enzymes, aligners, aggregation profiles, sample selection strategies, feature calling parameterizations
# Multiple versions for different -profile settings, containerization solutions

logfile="test_localPC.log"
logged="| tee -a $logfile"



log() {
    test_name="$1"
    command="$2"
    eval "echo \"---------------------\" $logged"
    eval "echo \"$test_name\" $logged"
    eval "echo \"$(date)\" $logged"
    eval "echo \"---------------------\" $logged"
    eval "echo \"$command\" $logged"
}

run() { 
    test_name="$1"
    command="$2"
    log "$test_name" "$command"
    eval "$command $logged" 
}

log "test_localPC" ""

# Define the initial part of the command in a variable
base_command="/usr/bin/time nextflow run contents/main.nf -profile localPC -c tests/nextflow.config -w testwork -resume -stub-run"

# Tests: Making bwa aligner index and other input resources
run "Run one HiC 3.0 (DdeI,DpnII) downsampled sample from Akgol 2021 with bwa." "$base_command --sampleFile vignettes/akgol2021/one_rep_bwa.tsv -params-file tests/onerep_bulk_hg38.yml"

# Tests: Tiered merge and deduplication
run "Run a tiered HiC 3.0 (DdeI,DpnII) downsampled condition from Akgol 2021 with bwa using resource files created in the last test." "$base_command --sampleFile vignettes/akgol2021/samples_bwa.tsv -params-file tests/tiered_bulk_hg38.yml"

# Tests: Using --fastqInterleaved to process interleaved fastq data, setting defaults via command line, BSBolt alignment, BSBolt aligner indexing
run "Run a bisulfite converted single-cell sample from snmC-seq3 from the Dekker lab with bsbolt, creating a new index." "$base_command -params-file tests/onerep_bulk_hg38.yml --fastqInterleaved \"vignettes/snmC-seq3/SRR13751346.fastq.gz\" --defaults.condition SRR13751346 --defaults.minMapq 10 --defaults.aligner bsbolt --defaults.assembly hg38 --defaults.restrictionEnzymes \"NlaIII,MboI\" --defaults.bwaFlags \"'-SP5Mp'\" --lastStep hicMatrix"



