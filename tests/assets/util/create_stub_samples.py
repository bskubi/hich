from pathlib import Path

import itertools as it
import polars as pl

cd = 10
br = 2
tr = 500
count = cd * br * tr
filename = f"{count}_stub_samples.tsv"
path = Path(filename).absolute()
samples = it.product(range(cd), range(br), range(tr))
df = pl.DataFrame(
    samples,
    schema = ["condition", "biorep", "techrep"],
    orient = 'row'
).with_columns(
    fastq1 = pl.lit("tests/assets/fastq/1k/1k_ERR1413593_1.fq.gz"),
    fastq2 = pl.lit("tests/assets/fastq/1k/1k_ERR1413593_2.fq.gz"),
    assembly = pl.lit("M129"),
    genomeReference = pl.lit("tests/assets/genomeReference/M129.fa.gz"),
    aligner = pl.lit("bwa-mem2"),
    restrictionEnzymes = pl.lit("HindIII")
)
print(path)
print(df)

df.write_csv(path, separator = "\t")


