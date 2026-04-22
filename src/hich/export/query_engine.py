from dataclasses import dataclass, field
import tiledb
from typing import Iterator, Literal
import pandas as pd
import duckdb
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq
import h5py
import click

from typing import Callable, List

def partition_data(
    df: pd.DataFrame, 
    partition_by: List[str], 
    strategy: Callable[[List[str], pd.DataFrame], None]
):
    """
    Isolates grouping logic. Generates path segments and dispatches 
    the subsetted DataFrame to the provided strategy callable.
    """
    cols_in_df = [c for c in partition_by if c in df.columns]
    data_cols = [c for c in df.columns if c not in cols_in_df]
    
    # Execute strategy immediately if no grouping is required
    if not cols_in_df:
        strategy(partition_by, df[data_cols])
        return

    for keys, subset in df.groupby(cols_in_df):
        key_vals = keys if isinstance(keys, tuple) else (keys,)
        key_map = dict(zip(cols_in_df, key_vals))
        
        path_segments = [
            str(key_map[item]) if item in key_map else str(item) 
            for item in partition_by
        ]
        
        strategy(path_segments, subset[data_cols])


class HDF5WriteStrategy:
    """
    Strategy class encapsulating HDF5 I/O and recarray conversion.
    """
    def __init__(self, file_obj: h5py.File, custom_dtypes: dict = None):
        self.file = file_obj
        self.custom_dtypes = custom_dtypes or {}

    def __call__(self, path_segments: List[str], df: pd.DataFrame):
        # 1. Convert to recarray
        if self.custom_dtypes:
            valid_dtypes = {k: v for k, v in self.custom_dtypes.items() if k in df.columns}
            if valid_dtypes:
                df = df.astype(valid_dtypes)
                
        payload = df.to_records(index=False)
        
        # 2. Traverse hierarchy
        node = self.file
        for segment in path_segments[:-1]:
            node = node.require_group(segment)
            
        dataset_name = path_segments[-1]
        
        # 3. Write or extend dataset
        if dataset_name in node:
            dset = node[dataset_name]
            curr_size = dset.shape[0]
            new_size = curr_size + len(payload)
            dset.resize((new_size,))
            dset[curr_size:] = payload
        else:
            node.create_dataset(
                dataset_name, 
                data=payload, 
                maxshape=(None,), 
                chunks=True
            )

class CSVWriteStrategy:
    """
    Strategy for writing chunked data to CSV/TSV format.
    Constructs a directory tree from path_segments.
    """
    def __init__(self, base_dir: str | Path, header: bool = True, sep: str = ",", ext: str | None = None):
        self.base_dir = Path(base_dir)
        self.header = header
        self.sep = sep
        if ext is None:
            self.ext = ".csv" if sep == "," else ".tsv"
        else:
            self.ext = ext

    def __call__(self, path_segments: List[str], df: pd.DataFrame) -> None:
        if not path_segments:
            file_path = self.base_dir.with_suffix(self.ext)
        else:
            file_path = self.base_dir.joinpath(*path_segments).with_suffix(self.ext)
            
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Determine write mode to support iterative chunk appending
        file_exists = file_path.exists()
        mode = "a" if file_exists else "w"
        write_header = self.header if not file_exists else False
        
        df.to_csv(
            file_path, 
            mode=mode, 
            sep=self.sep, 
            header=write_header, 
            index=False
        )


class ParquetWriteStrategy:
    """
    Strategy for writing chunked data to Parquet format.
    Requires 'fastparquet' engine for the append operations.
    """
    def __init__(self, base_dir: str | Path, ext: str = ".parquet"):
        self.base_dir = Path(base_dir)
        self.ext = ext

    def __call__(self, path_segments: List[str], df: pd.DataFrame) -> None:
        if not path_segments:
            file_path = self.base_dir.with_suffix(self.ext)
        else:
            file_path = self.base_dir.joinpath(*path_segments).with_suffix(self.ext)
            
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        append_mode = file_path.exists()
        
        df.to_parquet(
            file_path, 
            engine="fastparquet", 
            append=append_mode, 
            index=False
        )

import os
from pathlib import Path
import pandas as pd
import scipy.sparse as sp
from typing import List

class SciPySparseWriteStrategy:
    """
    Strategy for writing chunked data to SciPy CSR matrices stored as .npz files.
    """
    def __init__(self, base_dir: str | Path, custom_dtypes: dict = None, ext: str = ".npz"):
        self.base_dir = Path(base_dir)
        self.custom_dtypes = custom_dtypes or {}
        self.ext = ext

    def __call__(self, path_segments: List[str], df: pd.DataFrame) -> None:
        if self.custom_dtypes:
            valid_dtypes = {k: v for k, v in self.custom_dtypes.items() if k in df.columns}
            if valid_dtypes:
                df = df.astype(valid_dtypes)

        # Convert DataFrame to CSR
        chunk_csr = sp.csr_matrix(df.values)
        
        # Resolve target path
        if not path_segments:
            file_path = self.base_dir.with_suffix(".npz")
        else:
            file_path = self.base_dir.joinpath(*path_segments).with_suffix(".npz")
            
        file_path.parent.mkdir(parents=True, exist_ok=True)
        
        # Write or extend existing matrix
        if file_path.exists():
            existing_csr = sp.load_npz(file_path)
            combined_csr = sp.vstack([existing_csr, chunk_csr], format="csr")
            sp.save_npz(file_path, combined_csr)
        else:
            sp.save_npz(file_path, chunk_csr)

@dataclass
class ExportEngine:
    tiledb_array: tiledb.SparseArray
    query_names: list[str]
    query_cond: str | None = None
    query_order: Literal["global", "col-major", "row-major", "unordered"] = "col-major"

    @property
    def all_dim_names(self) -> list[str]:
        return [self.tiledb_array.dim(i).name for i in range(self.tiledb_array.ndim)]

    @property
    def all_attr_names(self) -> list[str]:
        return self.tiledb_array.attr_names
    
    @property
    def query_dim_names(self) -> list[str]:
        return [name for name in self.all_dim_names if name in self.query_names]
    
    @property
    def query_attr_names(self) -> list[str]:
        return [name for name in self.all_attr_names if name in self.query_names]
    
    @property
    def _query_order(self) -> str:
        match self.query_order:
            case "global": return  "G"
            case "col-major": return "F"
            case "row-major": return "C"
            case "unordered": return "U"
    
    @property
    def _query_names(self) -> list[str] | slice:
        return self.query_names or slice(None)

    @property
    def query(self) -> tiledb.Query:
        return self.tiledb_array.query(
            attrs=self.query_attr_names,
            cond=self.query_cond,
            dims=self.query_dim_names,
            order=self._query_order,
            return_incomplete=True
        )

    def query_tiledb(self) -> Iterator[pd.DataFrame]:
        for df_chunk in self.query.df[:]:
            yield df_chunk[self._query_names]

    def write_partition(self, strategy, partition_by: list[str]):
        for df_chunk in self.query.df[:]:
            partition_data(df_chunk[self._query_names], partition_by, strategy)
    
    def write_csv(self, target, header = False, sep = ",", **kwargs):
        for i, df in enumerate(self.query_tiledb()):
            write_header, mode = (header, "w") if i == 0 else (False, "a")
            
            df.to_csv(
                target, 
                sep=sep, 
                header=write_header, 
                mode=mode, 
                **kwargs
            )

    def write_parquet(self, target, **kwargs):
        writer = None
        for df in self.query_tiledb():
            table = pa.Table.from_pandas(df)
            if writer is None:
                writer = pq.ParquetWriter(target, table.schema)
            writer.write_table(table, **kwargs)
        
        if writer:
            writer.close()

@dataclass
class DuckDBExportEngine(ExportEngine):
    sql: list[str] = field(default_factory=list)
    duckdb_path: str = ":memory:"
    table_name: str = "QUERY_RESULT"
    vectors_per_chunk: int = 100_000

    def create_table_sql(self, df_name: str) -> str:
        return f"CREATE TABLE {self.table_name} AS SELECT * FROM {df_name}"
    
    def insert_into_table_sql(self, df_name: str) -> str:
        return f"INSERT INTO {self.table_name} SELECT * FROM {df_name}"        

    def query_tiledb(self) -> Iterator[pd.DataFrame]:
        with duckdb.connect(self.duckdb_path) as con:
            for i, df_chunk in enumerate(self.query.df[:]):
                df_chunk = df_chunk[self._query_names]
                if i == 0:
                    exec_result = con.execute(self.create_table_sql("df_chunk"))
                else:
                    exec_result = con.execute(self.insert_into_table_sql("df_chunk"))
            
            for sql in self.sql:
                exec_result = con.execute(sql)

            while True:
                # Fetch a single chunk
                df_chunk = exec_result.fetch_df_chunk(self.vectors_per_chunk)
                # DuckDB returns an empty DataFrame when the result set is exhausted
                if df_chunk.empty:
                    break
    
                yield df_chunk
            
    def write_partition(self, strategy, partition_by: list[str]):
        for df_chunk in self.query_tiledb():
            partition_data(df_chunk, partition_by, strategy)

    def write_csv(self, target, header = False, sep = ",", select: str | None = None, **kwargs):
        self.validate_con()

        select = select or self.table_name
        self.con.execute(f"COPY ({select}) TO {target} (FORMAT CSV, DELIMITER {sep}, header {str(header).lower()})")

    def write_parquet(self, target, select: str | None = None, **kwargs):
        self.validate_con()

        select = select or self.table_name
        self.con.execute(f"COPY ({select}) TO {target} (FORMAT PQRQUET)")
