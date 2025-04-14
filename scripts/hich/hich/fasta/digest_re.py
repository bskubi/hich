from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Restriction import RestrictionBatch
from pathlib import Path
import smart_open_with_pbgzip
from smart_open import smart_open
from typing import *
import itertools
import polars as pl

protocols = {
    'HiC 3.0': ['DpnII', 'DdeI'],
    'Hi-C 3.0': ['DpnII', 'DdeI'],
    'Arima Genome-Wide HiC+': ['DpnII', 'HinfI'],
    'Arima': ['DpnII', 'HinfI'],
    'Phase Proximo 2021+ Plant': ['DpnII'],
    'Phase Plant': ['DpnII'],
    'Phase Proximo 2021+ Animal': ['DpnII'],
    'Phase Animal': ['DpnII'],
    'Phase Proximo 2021+ Microbiome': ['Sau3AI', 'MluCI'],
    'Phase Microbiome': ['Sau3AI', 'MluCI'],
    'Phase Proximo 2021+ Human': ['DpnII'],
    'Phase Human': ['DpnII'],
    'Phase Proximo 2021+ Fungal': ['DpnII'],
    'Phase Fungal': ['DpnII']
}

def unique_enzymes(enzymes_and_protocols):
    """Extract unique entries from flattened iterable"""
    protocols_upper = {k.upper(): v for k, v in protocols.items()}
    enzymes = [protocols_upper.get(it.upper(), [it]) for it in enzymes_and_protocols]
    return list(set(itertools.chain.from_iterable(enzymes)))

def sequences(fasta: Path | str) -> Generator[SeqRecord, None, None]:
    """Yields sequences from fasta file (may be compressed)
    """
    handle = smart_open(fasta, "rt")
    yield from SeqIO.parse(handle, "fasta")

def pos_to_bed(chrom: str, pos: List[int], end: int, start: int = 1) -> pl.DataFrame:
    """Convert list of positions to a BED dataframe
    """
    # Ensure starts contains zero and not the last value
    # Ensure ends contains the last value and not zero
    starts = sorted(list(set([p for p in [start, *pos] if p != end])))
    ends = sorted(list(set([p for p in [*pos, start + end] if p != start]))) 
    chroms = [chrom]*len(starts)

    return (
        pl.DataFrame({
            "chrom": chroms,
            "start": starts,
            "end": ends
        })
    )


def digest_re(fasta: Path | str, enzymes_and_protocols: List[str]):
    enzymes = unique_enzymes(enzymes_and_protocols)
    try:
        restriction_batch = RestrictionBatch(enzymes)
    except Exception as e:
        print(f"Failed to create RestrictionBatch using enzymes list {enzymes}")
        print(e)

    for sequence in sequences(fasta):
        digest = restriction_batch.search(sequence.seq)
        pos = list(itertools.chain.from_iterable(digest.values()))
        yield pos_to_bed(sequence.id, pos, len(sequence.seq))

