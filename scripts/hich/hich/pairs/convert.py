"""Interconversions between pairs, pairsql, pairsparquet"""
from pathlib import Path
from typing import *
import os
from parse import parse

def rename(name, pattern1, pattern2):
    """Reformat name"""
    parsed = parse(pattern1, name)
    if parsed is None:
        raise ValueError(f"Could not parse path {name} using pattern {pattern1}")
    named = parsed.named
    return pattern2.format(**named)

def autodetect_extension(path: Path | str, default: Tuple[str, str] = (None, None)) -> Tuple[str, str]:
    """Detect unique format-specifying extension from path
    
    Returns:
        Tuple[str]: format and extension identified in path 
    
    Raises:
        Exception if multiple file format extensions are detected
    """
    if path is None:
        return None, None

    path = Path(path)
    suffixes = set(path.suffixes)

    autodetect = {
        "duckdb": [".duckdb", ".ddb", ".db", ".pairsduckdb"],
        "pairs": [".pairs"],
        "parquet": [".parquet", ".pq", ".pairsparquet"]
    }

    path_format = None
    format_extension = None
    for format, extensions in autodetect.items():
        detected_extension = set(extensions).intersection(suffixes)
        detected_extension = detected_extension.pop() if detected_extension else None

        repeated = path_format and detected_extension
        assert not repeated, f"Multiple file format extensions autodetected in {path}: {detected_extension} and {format_extension}."

        if detected_extension:
            path_format = format
            format_extension = detected_extension
    
    result = (path_format, format_extension) if path_format and format_extension else default

    assert result[0] and result[1], f"File format and extension not detected in {path} with suffixes {path.suffixes}. Options include {autodetect}. "

    return result


def walk_files(root_path: Path | str):
    """Walk through files in a directory"""
    for root, dirs, files in os.walk(root_path):
            if dirs:
                for dir in dirs:
                    yield from walk_files(root_path / dir)
            for file in files:
                yield Path(root) / file

def reorganize(out_path: Path | str, in_pattern: str = None, out_pattern: str = None, in_format: str = None, out_format: str = None, unlink: bool = True):
    """Walk through directory files and convert file formats and file names.
    
    Context: helper function for PairsSQL.partition_by
    """
    from hich.pairs.pairssql import PairsSQL
    assert not out_pattern or not Path(out_pattern).is_relative_to(out_path), f"--out-pattern {out_pattern} generates files in a subdirectory of OUT_PATH {out_path}, which is not allowed. Try using a different base directory in --out-pattern."
    out_path = Path(out_path)
    for partition_file in walk_files(out_path):
        if out_format != "parquet" or (in_pattern and out_pattern):
            if in_pattern and out_pattern:
                destination_file = Path(rename(str(partition_file), in_pattern, out_pattern))
            else:
                destination_file = Path(partition_file)

            destination_file.parent.mkdir(parents=True, exist_ok=True)

            if out_format == "parquet" or (not out_format and autodetect_extension(destination_file) == "parquet"):
                partition_file.rename(destination_file)
            else:
                db = PairsSQL.open(partition_file, out_format)

                db.write(destination_file, out_format)
                if partition_file.absolute() != destination_file.absolute() and unlink:
                    partition_file.unlink()