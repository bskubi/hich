mkdir -p bwa/mem bwa/mem2 bwameth/mem bwameth/mem2

bwa index -p hg38 hg38.fa.gz
mv hg38.ann hg38.bwt hg38.pac hg38.sa hg38.amb bwa/mem

bwa-mem2 index -p hg38 hg38.fa.gz
mv hg38.ann hg38.bwt hg38.pac hg38.sa hg38.amb bwa/mem2