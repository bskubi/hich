import tiledb as tdb
from hich.core import Category
import polars as pl
import numpy as np
import pysam

def read_seg(seg: pysam.AlignedSegment) -> pysam.AlignedSegment:
    return {
        "barcode": seg.get_tag("ZB") if seg.has_tag("ZB") else "",
        "readID": seg.query_name,
        "is_primary": not seg.is_secondary and not seg.is_supplementary,
        "is_mapped": seg.is_mapped,
        "is_duplicate": seg.is_duplicate
    }

def validate_observations(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df
        .with_columns(barcode = pl.col.barcode.str.replace("\.R$",""))
        .filter(pl.col.barcode != "", pl.col.context != "CN", ~pl.col.is_duplicate)
    )

def statistics() -> dict:
    return {
        "context": Category(
            lambda df: (
                df.select("barcode", met = pl.when(pl.col.mC).then(pl.lit("m")).otherwise(pl.lit("")), context = pl.col.context)
            ),
            pivot_on = ["met", "context"],
            pivot_index = ["barcode"],
            rename_pivot_on_pattern="{met}{context}"
        ),
        "reads": Category(
            lambda df: (
                (
                    df
                    .filter("is_primary")
                    .select("barcode", "readID", is_mapped = pl.when("is_mapped").then(pl.lit("Mapped")).otherwise(pl.lit("Unmapped")))
                    .unique()
                    .select("barcode", "is_mapped")
                )
            ),
            pivot_on = ["is_mapped"],
            pivot_index = ["barcode"]
        ),
    }

def observations_schema():
    return tdb.ArraySchema(
        domain = tdb.Domain(
            tdb.Dim(name="context", dtype="ascii", filters=fl_dict_zstd),
            tdb.Dim(name="chrom", dtype="ascii", filters=fl_dict_zstd),
            tdb.Dim(name="pos", domain=maximize_domain(pos_dt := np.int64, pos_tile := 10_000_000), dtype = pos_dt, tile = pos_tile, filters=fl_zstd),
            tdb.Dim(name="barcode", dtype="ascii", filters=fl_dict_zstd)
        ),
        attrs = [
            tdb.Attr(name="mC", dtype=np.bool_, filters=fl_zstd),
            tdb.Attr(name="strand", dtype="ascii", filters=fl_dict_zstd)
        ],
        sparse=True,
        allows_duplicates=True
    )

# Atomic compression methods
f_zstd = tdb.ZstdFilter(level=3)
f_dict = tdb.DictionaryFilter(level=3)

# Combined compression methods
fl_zstd = tdb.FilterList([f_zstd])
fl_dict_zstd = tdb.FilterList([f_dict, f_zstd])

# Maximize domain
maximize_domain = lambda dtype, tile: (np.iinfo(dtype).min, np.iinfo(dtype).max - tile - 1)
