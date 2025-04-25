import click
import smart_open_with_pbgzip
from smart_open import smart_open
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
@click.option("--memory-limit", type = int, default = None, help = "DuckDB memory limit in GB")
@click.option("--threads", type = int, default = None, help = "DuckDB thread limit")
@click.argument("in-path")
@click.argument("out-path")
@click.argument("partition_by", nargs=-1)
def partition(in_format, out_format, in_pattern, out_pattern, sql, squote, unlink, memory_limit, threads, in_path, out_path, partition_by):
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
        db.set_memory_limit(memory_limit)
        db.set_threads(threads)
        
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
@click.option("--out-path", default = None, help = "If supplied, changes are rewritten to this file, otherwise to stdout")
@click.option("--print-sql", default = False, is_flag = True, help = "Print SQL instead of running it")
@click.option("--memory-limit", type = str, default = None, help = "DuckDB memory limit in GB")
@click.option("--threads", type = int, default = None, help = "DuckDB thread limit")
@click.argument("sql")
@click.argument("in-path")
def sql(in_format, out_format, squote, out_path, print_sql, memory_limit, threads, sql, in_path):
    """Run a DuckDB SQL query on a .pairs-like file

    The 4DN .pairs format is ingested to '.pairssql' format using DuckDB, which has a `pairs` table having the same columns and names as the original .pairs file. Column types are autodetected through a full scan of the entire .pairs file. If outputting to .pairs, the header will be updated with any changed column names. If outputting to Parquet or DuckDB, the output will store the original .pairs header, either as a parquet kv metadata value "header" or the DuckDB table "metadata". The header will lack the #columns: line as this is generated on the fly when outputting to .pairs from the pairs table columns. 

    \b
    SQL: The DuckDB SQL query to run over file after ingesting to DuckDB. May also be a path to a file containing an SQL command.
    IN_PATH: Path to input file; use /dev/stdin to read from stdin

    \b
    Examples:
    \b
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
        # Load SQL from file
        sql_path = Path(sql)
        if sql_path.exists():
            sql = smart_open(sql_path).read()
    except:
        pass
    if squote:
        sql = sql.replace(squote, "'")
    if print_sql:
        print(sql)
    else:
        db = PairsSQL.open(in_path, in_format)
        db.set_memory_limit(memory_limit)
        db.set_threads(threads)
        try:
            if sql:
                db.conn.execute(sql)
        except Exception as e:
            print(f"SQL command failed on {in_path}:\n{sql}")
            print(e)
        db.write(out_path, out_format)

@pairs.command()
@click.option("--idx1", default = "rfrag1")
@click.option("--start1", default = "rfrag_start1")
@click.option("--end1", default = "rfrag_end1")
@click.option("--idx2", default = "rfrag2")
@click.option("--start2", default = "rfrag_start2")
@click.option("--end2", default = "rfrag_end2")
@click.option("--chrom1", default = "chrom1")
@click.option("--pos1", default = "pos1")
@click.option("--chrom2", default = "chrom2")
@click.option("--pos2", default = "pos2")
@click.option("--chunk-size", type = int, default = 1000000, show_default = True, help = "Number of lines to parse at a time")
@click.argument("bed-intervals")
@click.argument("in-path")
@click.argument("out-path")
def map_ends(idx1, start1, end1, idx2, start2, end2, chrom1, pos1, chrom2, pos2, chunk_size, bed_intervals, in_path, out_path):
    from hich.pairs.map_ends import map_ends as map_ends_HELPER
    from hich.pairs.header import PairsHeader
    import pandas as pd
    from io import StringIO
    header = PairsHeader(in_path)
    input_columns = header.columns.copy()
    header.columns += [idx1, start1, end1, idx2, start2, end2]
    intervals = pd.read_csv(
        bed_intervals, 
        delimiter = "\t", 
        names = ["chrom", "start", "end"], 
        dtype = {"chrom": str, "start": np.int64, "end": np.int64}
    )
    intervals["index"] = intervals.index
    intervals = pl.from_pandas(intervals)
    in_handle = smart_open(in_path, "rt")
    out_handle = smart_open(out_path, "wt")
    header.lines[-1] = f"#columns: " + " ".join(header.columns) + "\n"
    out_handle.write(header.text)
    pl.Config.set_tbl_cols(-1)

    for chunk in pd.read_csv(in_handle, delimiter="\t", skiprows=len(header.lines), names=input_columns, chunksize=chunk_size):
        chunk = pl.from_pandas(chunk)
        chunk = map_ends_HELPER(chunk, intervals, chrom1, pos1, chrom2, pos2, "chrom", "start", "end", idx1, start1, end1, idx2, start2, end2)
        buffer = StringIO()
        chunk.write_csv(buffer, include_header = False, separator = "\t")
        buffer.seek(0)
        out_handle.write(buffer.read())
