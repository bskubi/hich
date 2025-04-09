import click
import smart_open_with_pbgzip
import duckdb
import polars as pl
import shutil
from pathlib import Path
from hich.pairs.convert import walk_files, rename, reorganize
from hich.pairs.pairssql import PairsSQL
import os
from parse import parse
import numpy as np

@click.group()
def pairs():
    pass

@pairs.command()
@click.option("--in-format", type = click.Choice(["autodetect", "duckdb", "pairs", "parquet"]), default = "autodetect", help = "Input file format")
@click.option("--out-format", type = click.Choice(["autodetect", "duckdb", "pairs", "parquet"]), default = "autodetect", help = "Output file format")
@click.option("--in-pattern", type = str, default = None, help = "Python parse format for extracting names from original partitioned file.")
@click.option("--out-pattern", type = str, default = None, help = "python parse format for creating new filename. Output files should NOT be created in OUT_PATH directory.")
@click.option("--sql", type = str, default = None, help = "SQL to run on input file before partition.")
@click.option("--squote", default = "\"", help = "Replace this string with a single quote ' in the sql string")
@click.option("--unlink", is_flag = True, default = False, help = "Delete original partitioned file if renaming.")
@click.argument("in-path")
@click.argument("out-path")
@click.argument("partition_by", nargs=-1)
def partition(in_format, out_format, in_pattern, out_pattern, sql, squote, unlink, in_path, out_path, partition_by):
    """Split a .pairs-like file into multiple .pairs-like outputs

    \b
    IN_PATH: Path to input file to be partitioned (.pairs, .pairssql, .pairsparquet)
    OUT_PATH: Location where partitioned output files will be initially generated
    [PARTITION_BY]: Columns to partition the output; one output file generated per combination of values in these columns

    By default, all files are stored as a pairsparquet file named "data_0.parquet" in a partition-specific subdirectory of OUT_PATH. Subdirectories reflect a tree structure based on values of PARTITION_BY. The names of the first tier of subdirectories are values of the first column in PARTITION_BY, the second tier reflects values in the second column in PARTITION_BY, etc.

    \b
    Examples:
    Split to per-chromosome pairsparquet files in the directory structure output_dir/chrom1={chrom1_val}/chrom2={chrom2_val}/data_0.parquet:
        "hich pairs partition all_cells.pairs output_dir chrom1 chrom2"
    Convert outputs to .pairs format files named ./results/{chrom1_val}_{chrom2_val}.pairs:
        "hich pairs partition --in-pattern "output_dir/chrom1={chrom1}/chrom2={chrom2}/data_0.parquet" --out-pattern "results/{chrom1}_{chrom2}.pairs all_cells.pairs" output_dir chrom1 chrom2"
    Split by same vs. different chromosomes when that was not already labeled in the .pairs file:
        "hich pairs partition --sql "ALTER TABLE pairs ADD COLUMN same_chrom BOOLEAN; UPDATE pairs SET same_chrom = (chrom1 = chrom2)" all_cells.pairs output_dir same_chrom
    """
    if not partition_by:
        raise ValueError(f"No column names to partition {in_path} by were submitted.")

    try:
        db = PairsSQL().open(in_path, in_format)
    except Exception as e:
        print(f"Failed to open {in_path} with format {in_format}")
        raise e

    try:
        if sql:
            if squote:
                sql = sql.replace(squote, "'")
            db.conn.execute(sql)
    except Exception as e:
        print(f"Preliminary SQL query failed: {sql}")
        raise e

    try:
        db.partition_by(out_path, partition_by)
    except Exception as e:
        print(f"Failed to partition {in_path} by {partition_by} in output directory {out_path} ")
        raise e
    
    try:
        reorganize(out_path, in_pattern, out_pattern, in_format, out_format, unlink)
    except Exception as e:
        raise e

@pairs.command()
@click.option("--in-format", type = click.Choice(["autodetect", "duckdb", "pairs"]), default = "autodetect", help = "Input file format")
@click.option("--out-format", type = click.Choice(["autodetect", "duckdb", "parquet", "pairs", "tsv", "csv", "sql"]), default = "autodetect", help = "Output file format")
@click.option("--squote", default = "\"", help = "Replace this string with a single quote ' in the sql string")
@click.option("--out-path", default = None, help = "If supplied, changes are rewritten to a new file")
@click.option("--print-sql", default = False, is_flag = True, help = "Print SQL")
@click.argument("sql")
@click.argument("in-path")
def sql(in_format, out_format, squote, out_path, print_sql, sql, in_path):
    """Run a DuckDB SQL query on a .pairs-like file

    The 4DN .pairs format is ingested to '.pairssql' format using DuckDB, which has a `pairs` table having the same columns and names as the original .pairs file. Column types are autodetected through a full scan of the entire .pairs file. If outputting to .pairs, the header will be updated with any changed column names. If outputting to Parquet or DuckDB, the output will store the original .pairs header, either as a parquet kv metadata value "header" or the DuckDB table "metadata". The header will lack the #columns: line as this is generated on the fly when outputting to .pairs from the pairs table columns. 

    \b
    SQL: The DuckDB SQL query to run over file after ingesting to DuckDB. May also be a file containing an SQL command.
    IN_PATH: Path to input file
    OUT_PATH: .pairs-like filename where results are saved (updates IN_PATH_PAIRS if they are the same)

    \b
    Examples:

    Extract the substring of the readID column prior to the first ':' character and set as the value of the cellID column
        hich pairs sql "ALTER TABLE pairs ADD COLUMN cellID VARCHAR; UPDATE pairs SET cellID = regexp_extract(readID, \"(.*):(.*)\", 1);" no_cellID.pairs cellID_labeled.pairs
    Add a log10 distance strata with null values for transchromosomal or zero-distance contacts
        hich pairs sql "ALTER TABLE pairs ADD COLUMN distance INTEGER; UPDATE pairs SET distance = ROUND(LOG10(pos2 - pos1)) WHERE chrom1 == chrom2 AND pos1 != pos2;"
    Drop contacts mapping to different chromosomes
        hich pairs sql "DELETE FROM pairs WHERE chrom1 != chrom2;"
    Count number of contacts mapping to different distance strata:
        hich pairs sql "CREATE TEMP TABLE pairs_counts AS SELECT CAST(ROUND(LOG10(pos2-pos1)) AS INTEGER) A
S distance, COUNT(*) AS count FROM pairs WHERE pos1 != pos2 AND chrom1 == chrom2 GROUP BY distance; DROP TABLE pairs; CREATE TABLE pairs AS SELECT * FROM pairs_counts;"
    """
    try:
        sql_path = Path(sql)
        if sql_path.exists():
            sql = open(sql_path).read()
    except:
        pass
    if squote:
        sql = sql.replace(squote, "'")
    if print_sql:
        print(sql)
    db = PairsSQL.open(in_path, in_format)
    if sql:
        db.conn.execute(sql)
    if out_path:
        db.write(out_path, out_format)

@pairs.command()
@click.option("--in-format", type = click.Choice(["autodetect", "duckdb", "pairs"]), default = "autodetect", help = "Input file format")
@click.option("--out-format", type = click.Choice(["autodetect", "duckdb", "parquet", "pairs", "sql"]), default = "autodetect", help = "Output file format")
@click.option("--idx1", default = "idx1", show_default=True, help = "Label of first index")
@click.option("--start1", default = "start1", show_default=True, help = "Label of second BED interval")
@click.option("--end1", default = "end1", show_default=True, help = "Label of first BED interval")
@click.option("--idx2", default = "idx2", show_default=True, help = "Label of second BED index")
@click.option("--start2", default = "start2", show_default=True, help = "Label of second BED start")
@click.option("--end2", default = "end2", show_default=True, help = "Label of second BED end")
@click.option("--batch-size", default = 1000, show_default=True, help = "Number of blocks of 1024 rows at a time to treat as batch size")
@click.argument("in_path_pairs")
@click.argument("in_path_bed")
@click.argument("out_path")
def bin(in_format, out_format, idx1, start1, end1, idx2, start2, end2, batch_size, in_path_pairs, in_path_bed, out_path):
    """Bin contacts using a nonuniform partition

    IN_PATH_PAIRS: A .pairs-like file to compute intersections on the ends of its contacts
    IN_PATH_BED: A BED file defining the partition
    OUT_PATH: .pairs-like filename where results are saved (updates IN_PATH_PAIRS if they are the same)
    """   
    same_file = Path(in_path_pairs) == Path(out_path)
    assert not same_file or in_format == out_format, "Cannot change format while keeping same filename."

    in_db = PairsSQL.open(in_path_pairs, in_format)
    out_db = PairsSQL(path=":memory:")
    out_db.add_metadata(in_db.header)


    # Validate that bed file has no gaps or overlaps
    bed = duckdb.read_csv(in_path_bed, names=["chrom", "start", "end"]).pl().sort("chrom", "start")
    row1_iter = bed.iter_rows()
    row2_iter = bed.iter_rows()
    next(row2_iter)
    for row1, row2 in zip(row1_iter, row2_iter):
        chroms_match = row1[0] == row2[0]
        gap_or_overlap = row2[1] != row1[2]
        if chroms_match and gap_or_overlap:
            raise ValueError(f"Bed file {in_path_bed} does not define a partition. Gap at {row1} and {row2}.")
        

    def intersection_labels(pairs_chunk, bed, chrom, pos_col, idx_label, start_label, end_label):
        bed_chrom = bed.filter(pl.col("chrom") == chrom)

        if bed_chrom is None or len(bed_chrom) == 0:
            idx = [0]*len(pairs_chunk)
            start = [-1]*len(pairs_chunk)
            end = [-1]*len(pairs_chunk)
        elif not bed_chrom["start"][0] == 0:
            raise ValueError(f"For {chrom}, BED file is missing start position 0")
        elif pairs_chunk[pos_col].max() >= bed_chrom["end"][-1]:
            raise ValueError(f"For {chrom}, found a position {pairs_chunk[pos_col].max()} higher than max BED file partition ending at {bed_chrom['end'][-1]}")
        else:
            idx = np.searchsorted(bed_chrom["end"], pairs_chunk[pos_col], side = "right")
            start = bed_chrom["start"][idx]
            end = bed_chrom["end"][idx]

        return {idx_label: idx, start_label: start, end_label: end}

    for chunk in in_db.iter_chroms(batch_size):
        chrom1 = chunk["chrom1"][0]
        frag1 = intersection_labels(chunk, bed, chrom1, "pos1", idx1, start1, end1)
        
        chrom2 = chunk["chrom2"][0]
        frag2 = intersection_labels(chunk, bed, chrom2, "pos2", idx2, start2, end2)
        tagged_chunk = pl.concat([pl.from_pandas(chunk), pl.DataFrame(frag1), pl.DataFrame(frag2)], how = "horizontal")


        tagged_chunk_schema = tagged_chunk.filter(False)
        if Path(in_path_pairs) == Path(out_path):
            table = "pairs_temp"
        else:
            table = "pairs"
        out_db.conn.execute(
f"""
CREATE TABLE IF NOT EXISTS {table} AS
SELECT *
FROM tagged_chunk_schema;

INSERT INTO {table}
SELECT * FROM tagged_chunk;
"""
        )

    if same_file:
        out_db.conn.execute(
"""
DROP TABLE IF EXISTS pairs;

CREATE TABLE pairs AS SELECT * FROM pairs_temp;

DROP TABLE IF EXISTS pairs_temp;
"""
        )

    out_db.write(out_path, out_format)
