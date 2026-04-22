import tiledb as tdb
from hich.core import Category
import polars as pl
import numpy as np

def validate_observations(df: pl.DataFrame) -> pl.DataFrame:
    return (
        df.rename({"ZB1":"barcode"})
        .with_columns(
            barcode = pl.col.barcode.str.replace("\.R$",""),
            is_duplicate = ((pl.col.ZD1 > 0) | (pl.col.ZD2 > 0))
        )
        .filter(pl.col.barcode != "", ~pl.col.is_duplicate)
    )

def statistics():
    return {
        "mapq": Category(
            lambda df: (
                pl.concat([df.select("barcode", mapq = pl.col.mapq1), df.select("barcode", mapq = pl.col.mapq2)])
            ),
            pivot_on = "mapq"
        ),
        "pair_type": Category(lambda df: df.select("barcode", "pair_type"), pivot_on = "pair_type"),
        "cis_distance_strand": Category(
            lambda df: (
                df.filter(pl.col.pos1 != pl.col.pos2, pl.col.chrom1 == pl.col.chrom2, pl.col.chrom1 != "!")
                .select(
                    "barcode", 
                    distance = 2**((pl.col.pos2 - pl.col.pos1).log(2).round(0)).cast(pl.Int64), 
                    strand = pl.col.strand1 + pl.col.strand2
                )
            ),
            pivot_on = ["distance"],
            pivot_index = ["barcode", "strand"]
        ),
        "cis_trans": Category(
            lambda df: (
                df.filter(pl.col.chrom1 != "!", pl.col.chrom2 != "!")
                .select("barcode", cis_trans = pl.when(pl.col.chrom1 == pl.col.chrom2).then(pl.lit("Cis")).otherwise(pl.lit("Trans")))
            ),
            pivot_on = "cis_trans",
            pivot_index = "barcode"
        )
    }

def observations_schema():
    return tdb.ArraySchema(
        domain = tdb.Domain(
            tdb.Dim(name="chrom1", dtype="ascii", filters=fl_dict_zstd),
            tdb.Dim(name="chrom2", dtype="ascii", filters=fl_dict_zstd),
            tdb.Dim(name="pos1", domain=maximize_domain(pos1_dt := np.int64, pos_tile := 10_000_000), dtype = pos1_dt, tile = pos_tile, filters=fl_zstd),
            tdb.Dim(name="pos2", domain=maximize_domain(pos2_dt := np.int64, pos_tile := 10_000_000), dtype = pos2_dt, tile = pos_tile, filters=fl_zstd),
            tdb.Dim(name="barcode", dtype="ascii", filters=fl_dict_zstd)
        ),
        attrs = [
            tdb.Attr(name="strand1", dtype="ascii", filters=fl_dict_zstd),
            tdb.Attr(name="strand2", dtype="ascii", filters=fl_dict_zstd)
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
