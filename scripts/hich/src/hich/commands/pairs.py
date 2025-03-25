import click
import duckdb as ddb
import smart_open_with_pbgzip
import duckdb
import shutil
from pathlib import Path
from hich.pairs.convert import walk_files, rename, reorganize
from hich.pairs.pairssql import PairsSQL
import os
from parse import parse

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
    """Define a partition to split a .pairs-like file into multiple outputs

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
    db = PairsSQL().open(in_path, in_format)
    if sql:
        if squote:
            sql = sql.replace(squote, "'")
        db.conn.execute(sql)
    
    db.partition_by(out_path, partition_by)
    
    reorganize(out_path, in_pattern, out_pattern, in_format, out_format, unlink)

@pairs.command()
@click.option("--in-format", type = click.Choice(["autodetect", "duckdb", "pairs"]), default = "autodetect", help = "Input file format")
@click.option("--out-format", type = click.Choice(["autodetect", "duckdb", "parquet", "pairs", "sql"]), default = "autodetect", help = "Output file format")
@click.option("--squote", default = "\"", help = "Replace this string with a single quote ' in the sql string")
@click.argument("sql")

@click.argument("in-path")
@click.argument("out-path")
def sql(in_format, out_format, squote, sql, in_path, out_path):
    """Run a DuckDB SQL query on a 4DN .pairs file or Parquet or DuckDB '.pairssql' file.

    The 4DN .pairs format is ingested to '.pairssql' format using DuckDB, which has a `pairs` table having the same columns and names as the original .pairs file. Column types are autodetected through a full scan of the entire .pairs file. If outputting to .pairs, the header will be updated with any changed column names. If outputting to Parquet or DuckDB, the output will store the original .pairs header, either as a parquet kv metadata value "header" or the DuckDB table "metadata". The header will lack the #columns: line as this is generated on the fly when outputting to .pairs from the pairs table columns. 

    \b
    SQL: The DuckDB SQL query to run over file after ingesting to DuckDB
    IN_PATH: Path to input file
    OUT_PATH: Path to output file where the results are saved

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
    if squote:
        sql = sql.replace(squote, "'")
    db = PairsSQL.open(in_path, in_format)
    if sql:
        db.conn.execute(sql)
    db.write(out_path, out_format)
