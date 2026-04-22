mkdir -p bwa/mem bwa/mem2 bwameth/mem bwameth/mem2

bwa mem index -p hg38 hg38.fa.gz
for f in hg38.fa.gz*; do
    mv "$f" "${f/hg38.fa.gz/hg38}"
done
mv hg38.ann hg38.bwt hg38.pac hg38.sa hg38.amb bwa

bwa-mem2 index -p hg38 hg38.fa.gz
for f in hg38.fa.gz*; do
    mv "$f" "${f/hg38.fa.gz/hg38}"
done
mv hg38.ann hg38.bwt hg38.pac hg38.sa hg38.amb bwa