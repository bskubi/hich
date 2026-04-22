from hich.core.statistics import *

import polars as pl
from hich import Category, Engine, Result, SAM_Engine
from typing import Iterator
import tiledb as tdb
import pysam
import time
import subprocess
import io
from pathlib import Path

import polars as pl
import subprocess
from typing import Any
from dataclasses import dataclass
import click

@dataclass
class Extract3DPairsEngine(SAM_Engine):
    command: str | list[str] | None = None

    @classmethod
    def cli_command(cls, f):
        decorators = [
            click.option("--stats-dir", default="./"),
            click.option("--n-workers", default=1),
            click.option("--worker-n-procs", default=1),
            click.option("--batch-size-bytes", default=32*1024*1024),
            click.argument("config"),
            click.argument("sam"),
            click.argument("observations"),
            click.argument("command")
        ]
        for d in decorators[::-1]:
            f = d(f)
        return f

    def SAM_batches(self, sam_stream: dict) -> dict:
        return {"command": self.command, **sam_stream}

    @classmethod
    def SAM_worker(cls, sam_stream: bytes, command: str | list[str]) -> Any:
        start_time = time.time()
        process = subprocess.run(
            command,
            shell = isinstance(command, str),
            input=sam_stream,
            capture_output=True,
            check=True,
            text=False
        )
        
        output_bytes = process.stdout
        if not output_bytes:
            return Result(pl.DataFrame())

        columns = None
        col_marker = b"#columns:"
        start_idx = output_bytes.find(col_marker)
        
        if start_idx != -1:
            end_idx = output_bytes.find(b"\n", start_idx)
            header_line = output_bytes[start_idx:end_idx].decode('ascii')
            columns = header_line.split()[1:]

        df = pl.read_csv(
            output_bytes,
            separator="\t",
            comment_prefix="#",
            has_header=False,
            new_columns=columns,
            n_threads=1 
        )[1:-1]
        

        return Result(df)
    

@click.command
@Extract3DPairsEngine.cli_command
def extract_3d_pairs(stats_dir, n_workers, worker_n_procs, batch_size_bytes, config, sam, observations, command):
    engine = Extract3DPairsEngine(
        statistics_dir=Path(stats_dir), 
        command = command, 
        max_workers=n_workers,
        samtools_view_cpus=worker_n_procs, 
        batch_size_bytes = batch_size_bytes
    )
    engine.load_config_py(config)
    engine.run(sam, observations, False)