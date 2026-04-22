from hich.core.statistics import *

import polars as pl
from hich import Category, Engine, Result, SAM_Engine
from typing import Iterator
import tiledb as tdb
import pysam
import time
import subprocess
import io
from pydantic import BaseModel
import polars as pl
import subprocess
from typing import Any, Self, Callable
from hich.pysam_utils import AlignmentFile
from collections import defaultdict
import numpy as np
import pyarrow as pa
import pysam
from pathlib import Path

from dataclasses import dataclass, field
from metex.metex import AlignedSegmentMethylation, extract_methylation_bwameth

from shared_ref.shared_ref import SharedRef
import click

class ExtractMethylationParams(BaseModel):
    # Offset of context with respect to cytosine (negative = 3', positive = 5')
    context_offset: int = 0

    # Length of context starting from pos + context_offset in 5' direction 
    context_length: int = 2

    # If context extends off chromosome, pad it with this character
    context_pad: str = '\x00'

    # Map exact contexts (values) to output contexts (keys)
    # Exact contexts not given are not edited. 
    context_names: dict[str, list[str]] = {
            "CH": ["CA", "CC", "CT"],
            "CHH": ["CAA", "CAC", "CAT"
                    "CCA", "CCC", "CCT",
                    "CTA", "CTC", "CTT"],
            "CHG": ["CAG", "CCG", "CCT"]
        }

    # Which contexts (after mapping with context_names) to filter out
    # for excessively high methylation 
    high_mC_context: list[str] = ["CH"]

    # Minimum number of observed cytosines to consider filtering out
    # for high mC
    high_mC_count: int = 4

    # If passing high_mC_count, minimum fraction of mC on read to filter it out.
    high_mC_frac: float = .5

class ExtractMethylationReference(BaseModel, arbitrary_types_allowed=True):
    """
    Parameters for orchestrating DNA methylation extraction
    """
    ref_path: str
    limit_chroms: list[str] = None
    ignore_chroms: list[str] = None
    ref_forward: SharedRef | int | None = None
    ref_reverse: SharedRef | int | None = None

    def load_ref(self) -> Self:
        """
        Load ref_forward and ref_reverse from fasta
        """
        # Setup variables
        limit_chroms = self.limit_chroms
        ignore_chroms = self.ignore_chroms
        ref_fa = pysam.FastxFile(self.ref_path)
        ref_forward_dict = {}
        ref_reverse_dict = {}
        complement = str.maketrans("ACTGN", "TGACN")
        found = {ch: False for ch in limit_chroms}

        for contig in ref_fa:
            # Ignore chroms (mainly for low-RAM debugging)
            if limit_chroms and contig.name not in limit_chroms:
                continue
            if ignore_chroms and contig.name in ignore_chroms:
                continue

            # Store {config_name: contig_seq} for forward, reference strand
            ref_forward_dict[contig.name] = contig.sequence.upper()
            ref_reverse_dict[contig.name] = contig.sequence.upper().translate(complement)
            
            # Halt early when all limit chroms found
            if contig.name in found:
                found[contig.name] = True
            if all(found.values()):
                break
        
        # SharedRef stores reference genome in shared memory
        # across all parallel processes.
        self.ref_forward = SharedRef(ref=ref_forward_dict)
        self.ref_reverse = SharedRef(ref=ref_reverse_dict)

        return self

    def dehydrate(self) -> Self:
        """Convert objects to IDs for serialization."""
        return self.model_copy(
            update={
                "ref_forward": self.ref_forward.id if hasattr(self.ref_forward, 'id') else self.ref_forward,
                "ref_reverse": self.ref_reverse.id if hasattr(self.ref_reverse, 'id') else self.ref_reverse
            }
        )
    
    def hydrate(self) -> Self:
        """Convert IDs back to objects."""
        return self.model_copy(
            update={
                "ref_forward": SharedRef(id=self.ref_forward) if isinstance(self.ref_forward, dict) else self.ref_forward,
                "ref_reverse": SharedRef(id=self.ref_reverse) if isinstance(self.ref_reverse, dict) else self.ref_reverse
            }
        )

@dataclass
class ExtractDNAMethylationEngine(SAM_Engine):
    reference: ExtractMethylationReference | None = None
    params: ExtractMethylationParams = field(default_factory=lambda: ExtractMethylationParams())
    read_seg: Callable[[pysam.AlignedSegment], pysam.AlignedSegment] | None = None

    @classmethod
    def cli_command(cls, f):
        decorators = [
            click.option("--stats-dir", default="./"),
            click.option("--limit-chroms", "--lc", multiple=True),
            click.option("--ignore-chroms", "--ic", multiple=True),
            click.option("--n-workers", default=1),
            click.option("--worker-n-procs", default=1),
            click.option("--batch-size-bytes", default=32*1024*1024),
            click.argument("config"),
            click.argument("sam"),
            click.argument("ref"),
            click.argument("observations")
        ]
        for d in decorators[::-1]:
            f = d(f)
        return f

    def SAM_batches(self, sam_stream: dict) -> dict:
        try: self.read_seg = self.ns.read_seg
        except: ...
        
        batch = dict(
            **sam_stream,
            params = self.params,
            reference = self.reference.dehydrate(),
            read_seg = self.read_seg
        )

        return batch
    
    @classmethod
    def SAM_worker(cls, sam_stream: bytes, params: ExtractMethylationParams, reference: ExtractMethylationReference, read_seg: Callable[[pysam.AlignedSegment], pysam.AlignedSegment]) -> Result:
        """
        Call methylation from 
        Arguments:
            sam_stream (str): SAM-format (plaintext) Header + reads for batch 
            plan (MetPlan): Configures how methylation will be called
        """
        # Connect to shared-memory reference genomes
        reference = reference.hydrate()

        # Create AlignmentFile from SAM batch
        file = AlignmentFile.from_text(sam_stream.decode())

        # Metex and user data are accumulated as list of per-seg positional arrays
        # User data has precedence over metex col names
        met_cols = ["context", "strand", "pos", "mC"]
        metex_data: defaultdict[str, list[np.ndarray]] = defaultdict(list)
        extra_data: defaultdict[str, list[np.ndarray]] = defaultdict(list)
        chrom_names = set(reference.ref_forward.contig_names + reference.ref_reverse.contig_names)

        # Extract params from user config or get default
        ctx_offset = params.context_offset
        ctx_length = params.context_length
        ctx_pad = params.context_pad
        ctx_rename = params.context_names
        high_mC_context = params.high_mC_context
        high_mC_count = params.high_mC_count
        high_mC_frac = params.high_mC_frac

        for seg in file:
            chrom = seg.reference_name
            if not chrom in chrom_names:
                continue

            aligned_pairs = seg.get_aligned_pairs(matches_only=True)
            if not aligned_pairs:
                continue

            met = extract_methylation_bwameth(
                seq=seg.query_sequence,
                is_read1=seg.is_read1,
                is_forward=seg.is_forward,
                aligned_pairs=aligned_pairs,
                ref_forward=reference.ref_forward[chrom],
                ref_reverse=reference.ref_reverse[chrom],
                ctx_offset=ctx_offset,
                ctx_length=ctx_length,
                ctx_pad=ctx_pad
            )
            
            count = met.count
            
            # Rename contexts
            for new_name, old_names in ctx_rename.items():
                met.context[np.isin(met.context, old_names)] = new_name

            # Add context, pos, strand, mC
            [metex_data[col].append(getattr(met, col)) for col in met_cols if col not in extra_data]
            
            # Add chrom name
            if "chrom" not in extra_data:
                metex_data["chrom"].append(np.repeat(chrom, count))

            # Add high mC filter
            if "high_mC" not in extra_data:
                high_mC_count = np.isin(met.context, high_mC_context).sum()
                high_mC = (
                    high_mC_count > high_mC_count 
                    and high_mC_count / count > high_mC_frac
                ) if count else False
                metex_data["high_mC"].append(np.repeat(high_mC, count))

            # Add user-defined extra cols
            if isinstance(read_seg, Callable):
                seg_extra_data = read_seg(seg)
                [extra_data[col].append(np.repeat(val, count)) for col, val in seg_extra_data.items()]

        # Create combined dataframe, with extra_data taking precedence
        # Note: don't use ChainMap here, it's incompatible with defaultdict
        merged_data = dict(metex_data) | dict(extra_data)
        arrow_dict = {}

        for col, arrays in merged_data.items():
            # pa.array is zero-copy for C-contiguous numpy arrays of native dtypes,
            # but a rechunk is forced if the DF is large enough.
            arrow_arrays = [pa.array(arr) for arr in arrays]
            arrow_dict[col] = pa.chunked_array(arrow_arrays)

        return Result(pl.from_arrow(pa.Table.from_pydict(arrow_dict)))
    
@click.command
@ExtractDNAMethylationEngine.cli_command
def extract_dna_methylation(stats_dir, limit_chroms, ignore_chroms, n_workers, worker_n_procs, batch_size_bytes, config, sam, ref, observations):
    reference = ExtractMethylationReference(
        ref_path = ref,
        limit_chroms=limit_chroms,
        ignore_chroms=ignore_chroms
    ).load_ref()
    engine = ExtractDNAMethylationEngine(
        statistics_dir=Path(stats_dir),
        max_workers=n_workers,
        reference = reference, 
        samtools_view_cpus=worker_n_procs, 
        batch_size_bytes = batch_size_bytes
    )
    engine.load_config_py(config)
    engine.run(sam, observations, False)