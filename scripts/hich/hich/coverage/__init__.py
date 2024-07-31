from itertools import combinations, product
from hich.parse.pairs_parser import PairsParser
import polars as pl
from dataclasses import dataclass, field
from typing import Tuple, List, Callable
from more_itertools import collapse
from numpy.random import default_rng
from hich.coverage.samheader_coverage import SamheaderCoverage

"""
Compute strata on each pairs file
Join on chrom1 chrom2 stratum

Strategies:
    - To given strata fractions
    - To mean strata fractions, min column total
    - To min column total, random sample

Create a way to randomly sample a specific number of pairs
    - Computing the strata gives us the total number to sample for each file
    - Use interpolation sampling since we know the total number of reads of each
      stratum in the file
"""

def drop_unmapped(df):
    return df.filter((pl.col("chrom1") != "!") & (pl.col("chrom2") != "!"))

def add_cis(df):
    return df.with_columns((pl.col("chrom1") == pl.col("chrom2")).alias("cis"))

def add_distance(df):
    return df.with_columns(
           pl.when(pl.col("cis"))
             .then((pl.col("pos2") - pl.col("pos1")).abs())
             .otherwise(-2)
             .alias("distance"))

def add_stratum(df, strata):
    strata_sorted = pl.Series([-1] + strata).sort()
    strata_labels = pl.Series([-1] + strata + [strata[-1] + 1])
    strata_indices = strata_sorted.search_sorted(df["distance"])
    labels = strata_labels[strata_indices]
    return df.with_columns(labels.alias("stratum"))

def add_to_count(total_counts, df, file):
    new_counts = df.group_by(["chrom1", "chrom2", "stratum"]).count()
    new_counts = new_counts.rename({"count":file})
    if total_counts is None:
        return new_counts
    return total_counts.join(new_counts,
                             on = ["chrom1", "chrom2", "stratum"],
                             how = "outer",
                             coalesce = True) \
                       .fill_null(0) \
                       .with_columns(pl.col(file) + pl.col(f"{file}_right")) \
                       .drop(f"{file}_right")

def with_stratum(df, strata):
    df = drop_unmapped(df)
    df = add_cis(df)
    df = add_distance(df)
    return add_stratum(df, strata)

def strata(file, strata):
    """Compute hic strata for a .pairs file; strata 0 is trans"""
    parser = PairsParser(file)
    

    count = None
    for df in parser.batch_iter(10000):
        df = with_stratum(df, strata)
        count = add_to_count(count, df, file)
    return count

def combine_strata(dfs):
    combined = dfs[0]
    for df in dfs[1:]:
        combined = combined.join(df,
                                 on = ["chrom1", "chrom2", "stratum"],
                                 how = "outer",
                                 coalesce = True)
    combined = combined.fill_null(0)
    return combined

@dataclass
class StratumSampler:
    N_records: int
    n_to_sample: int
    t_viewed: int = 0
    m_sampled: int = 0
    U_random: list[float] = field(default_factory=list)
    samples: list[bool] = field(default_factory=list)

    def batch(self):
        """Add a number of sample values equal to the length of self.U_random"""
        for u in self.U_random:
            if (self.N_records - self.t_viewed) * u < self.n_to_sample - self.m_sampled:
                self.m_sampled += 1
                self.samples.append(True)
            else:
                self.samples.append(False)
            self.t_viewed += 1
    
    def sample(self, n):
        """Pop and return up to the first n items from self.samples"""
        popped = self.samples[:n]
        self.samples = self.samples[n:]
        return popped


class MultiSampler:
    def __init__(self, df, id_cols, N_records_col, n_to_sample_col, seed = 42):
        self.setup(df, id_cols, N_records_col, n_to_sample_col, seed)

    def setup(self, df, id_cols, N_records_col, n_to_sample_col, seed = 42):
        sel = collapse([id_cols, N_records_col, n_to_sample_col])
        df = df.select(sel)

        N_col = len(id_cols)
        n_col = N_col + 1

        self.id_cols = id_cols
        self.rng = default_rng(seed=seed)
        self.samplers = {}
        for row in df.iter_rows():
            sampler_id = row[:len(id_cols)]
            N = row[N_col]
            n = row[n_col]
            self.samplers[sampler_id] = StratumSampler(N, n)
    
    def sample(self, df):
        groups = df.with_row_index() \
                   .group_by(self.id_cols)
        sampled = []
        for name, data in groups:
            n = len(data)
            sample = self.batch(name, n)
            sampled.append(data.filter(sample))
        
        index_col = df.columns[0]
        sampled = pl.concat(sampled).sort(by=index_col)
        return sampled
    
    def batch(self, sampler_id, n):
        self.samplers[sampler_id].U_random += self.rng.random(n).tolist()
        self.samplers[sampler_id].batch()
        return self.samplers[sampler_id].sample(n)


def selection_sample(in_file, out_file, strata, to_sample, sampler_id_cols, N_col, n_col, batch_size = 10000, seed = 42):
    sampler = MultiSampler(to_sample, sampler_id_cols, N_col, n_col)
    parser = PairsParser(in_file)
    for df in parser.batch_iter(batch_size):
        orig_cols = df.columns
        df = with_stratum(df, strata)
        sample = sampler.sample(df).select(orig_cols)

        parser.write_append(out_file, sample, header_end = SamheaderCoverage())

