mkdir -p assets/downloads/mcool/full
mkdir -p assets/downloads/genomeReference
mkdir -p assets/downloads/fastq

wget --no-clobber -O "assets/downloads/mcool/full/GSM4604274_629.iced.mcool" "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4604nnn/GSM4604274/suppl/GSM4604274%5F629.iced.mcool"
wget --no-clobber -O "assets/downloads/mcool/full/GSM4604271_576.iced.mcool" "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4604nnn/GSM4604271/suppl/GSM4604271%5F576.iced.mcool"
wget --no-clobber -O "assets/downloads/mcool/full/GSM4604272_103.iced.mcool" "https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM4604nnn/GSM4604272/suppl/GSM4604272%5F103.iced.mcool"
wget --no-clobber -O "assets/downloads/genomeReference/hg38.fa.gz" "https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
wget --no-clobber -O "assets/downloads/fastq/U54-HFFc6-DSG-DdeI-DpnII-20190711-R1-T1_S12_L004_R1_001.fastq.gz" "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/1889f325-776f-4dc3-81b2-7c22e19ad516/4DNFI7G518XA.fastq.gz"
wget --no-clobber -O "assets/downloads/fastq/U54-HFFc6-DSG-DdeI-DpnII-20190711-R1-T1_S12_L004_R2_001.fastq.gz" "https://4dn-open-data-public.s3.amazonaws.com/fourfront-webprod/files/4255b8f8-2581-4a84-a25a-705cafa8a2b2/4DNFI9H22RJ3.fastq.gz"
wget --no-clobber -O "assets/downloads/genomeReference/hg38.fa.gz" "https://www.encodeproject.org/files/GRCh38_no_alt_analysis_set_GCA_000001405.15/@@download/GRCh38_no_alt_analysis_set_GCA_000001405.15.fasta.gz"
wget --no-clobber -O "assets/downloads/genomeReference/mm10.fa.gz" "https://www.encodeproject.org/files/mm10_no_alt_analysis_set_ENCODE/@@download/mm10_no_alt_analysis_set_ENCODE.fasta.gz"

mkdir -p assets/index/bwa
cd assets/index/bwa
ln -s ../../downloads/genomeReference/hg38.fa.gz hg38
ln -s ../../downloads/genomeReference/mm10.fa.gz mm10
bwa index hg38
bwa index mm10
cd ../../..

mkdir -p assets/index/bwa-mem2
cd assets/index/bwa-mem2
ln -s ../../downloads/genomeReference/hg38.fa.gz hg38
ln -s ../../downloads/genomeReference/mm10.fa.gz mm10
bwa-mem2 index hg38
bwa-mem2 index mm10
cd ../../..

mkdir -p assets/index/bwameth
cd assets/index/bwameth
ln -s ../../downloads/genomeReference/mm10.fa.gz mm10
bwameth.py index mm10
cd ../../..

mkdir -p assets/index/bwameth-mem2
cd assets/index/bwameth-mem2
ln -s ../../downloads/genomeReference/mm10.fa.gz mm10
bwameth.py index-mem2 mm10
cd ../../..

mkdir -p assets/fragmentIndex
hich fasta re-fragments assets/downloads/genomeReference/hg38.fa.gz assets/fragmentIndex/hg38_DpnII,DdeI.bed DpnII DdeI
hich fasta re-fragments assets/downloads/genomeReference/mm10.fa.gz assets/fragmentIndex/mm10_DpnII,DdeI.bed DpnII DdeI