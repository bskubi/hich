from itertools import combinations, product
from hich.parse.pairs_parser import PairsParser
import polars as pl

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


def strata(file, strata):
    """Compute hic strata for a .pairs file; strata 0 is trans"""
    parser = PairsParser(file)
    

    count = None
    for df in parser.batch_iter(10000):
        df = drop_unmapped(df)
        df = add_cis(df)
        df = add_distance(df)
        df = add_stratum(df, strata)
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


