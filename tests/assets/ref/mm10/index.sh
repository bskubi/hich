python -m bwameth index mm10.fa.gz
for f in mm10.fa.gz*; do
    mv "$f" "${f/mm10.fa.gz/mm10}"
done
mv mm10.bwameth* bwameth/mem

python -m bwameth index-mem2 mm10.fa.gz
for f in mm10.fa.gz*; do
    mv "$f" "${f/mm10.fa.gz/mm10}"
done
mv mm10.bwameth* bwameth/mem2