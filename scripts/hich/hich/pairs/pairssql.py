""".pairssql is a DuckDB version of 4DN .pairs format"""
import smart_open_with_pbgzip
from smart_open import smart_open
import duckdb
from hich.pairs.convert import *
from hich.pairs.header import * 
from dataclasses import dataclass
import sys

@dataclass
class PairsSQL:
    path: Path = None
    _conn: duckdb.DuckDBPyConnection = None
    memory_limit: str = None

    def set_memory_limit(self, memory_limit: str | int | float):
        try:
            memory_limit = int(memory_limit)
        except:
            pass
        if type(memory_limit) in [float, int]:
            memory_limit = f"{memory_limit}GB"
        elif type(memory_limit) == str:
            memory_limit = "".join(memory_limit.split())

        if memory_limit is not None:
            self.conn.execute("SET memory_limit = $memory_limit", {"memory_limit": memory_limit})

    def set_threads(self, threads: int | None):
        if threads is not None:
            self.conn.execute("SET threads = $threads", {"threads": threads})

    def from_pairs(self, in_path: Path, out_path: Path | str = ":memory:"):
        self.connect(out_path)
        
        header = PairsHeader(in_path)
        
        
        self.path = out_path
        sql = (
f"""
DROP TABLE IF EXISTS pairs;
DROP TABLE IF EXISTS metadata;

CREATE TABLE pairs AS
SELECT * FROM read_csv($in_path, names={header.columns}, delim='\t', skip={len(header.lines)}, sample_size=-1, header=false);
"""
        )

        self.conn.execute(
        sql,
        {"in_path": str(in_path)}
        )

        self.add_metadata(header)
        return self
    
    def from_duckdb(self, in_path: Path | str, out_path: Path | str = ":memory:"):
        self.connect(out_path)
        if in_path != out_path:
            self.conn.execute(
f"""
ATTACH '{in_path}' AS db1 (READ_ONLY);
ATTACH '{out_path}' AS db2;
DROP TABLE IF EXISTS pairs;
DROP TABLE IF EXISTS metadata;
COPY FROM DATABASE db1 TO db2;
DETACH db1;
USE db2;
"""
            )
        return self

    def write_duckdb(self, out_path: Path | str):
        self.conn.execute(
f"""
ATTACH '{out_path}' AS write_duckdb;
DROP TABLE IF EXISTS write_duckdb.pairs;
DROP TABLE IF EXISTS write_duckdb.metadata;
COPY FROM DATABASE {self.current_database} TO write_duckdb;
DETACH write_duckdb;
"""
        )
    
    def from_parquet(self, in_path: Path | str, out_path: Path | str = ":memory:"):
        self.connect(out_path)
        self.conn.execute(
f"""
CREATE TABLE pairs AS SELECT * FROM '{in_path}';
"""
        )
        
        metadata = self.conn.execute(f"SELECT * FROM parquet_kv_metadata('{in_path}')").fetchall()
        header_text = [row[2].decode('utf-8') for row in metadata if row[1].decode('utf-8') == "header"][0]
        header = PairsHeader(header_text, from_text = True, columns = self.columns, replace_columns = True)
        self.add_metadata(header)
        return self

    def add_metadata(self, header: PairsHeader):
        self.conn.execute(
f"""
CREATE TABLE IF NOT EXISTS metadata (
    table_name VARCHAR,
    header VARCHAR
);

DELETE FROM metadata WHERE table_name = 'pairs';

INSERT INTO metadata VALUES ('pairs', '{header.non_columns_text}');
"""
        )

    def write_pairs(self, out_path: Path | str | None, vector_multiple = 1000):
        handle = smart_open(out_path, "w") if out_path is not None else sys.stdout

        handle.write(self.header.text)
        query = self.conn.execute("SELECT * FROM pairs")

        if vector_multiple is None:
            query.df().to_csv(handle, sep="\t", header=False, index=False, compression=None)
        else:
            while (chunk := query.fetch_df_chunk(vector_multiple)) is not None:
                if chunk.empty:
                    break
                else:
                    chunk.to_csv(handle, sep="\t", header=False, index=False, compression=None)
    

    def write_parquet(self, out_path: Path | str):
        self.conn.execute(
f"""
COPY (SELECT * FROM pairs)
TO '{out_path}'
(
    FORMAT PARQUET,
    KV_METADATA {{
        header: '{self.header.text}'
    }}
);
"""
        )

    def write(
            self,
            out_path: Path | str,
            out_format: str = "autodetect"
    ) -> duckdb.DuckDBPyConnection:
        """Write to pairs, duckdb, or parquet"""

        out_path = Path(out_path) if out_path else None

        if not out_format or out_format == "autodetect":
            out_format, _ = autodetect_extension(out_path)

        if out_format == "duckdb":
            self.write_duckdb(out_path)
        elif out_format == "parquet":
            self.write_parquet(out_path)
        elif out_format == "pairs" or out_path is None:
            self.write_pairs(out_path)
        else:
            raise Exception(f"Format {out_format} not recognized. Should be one of duckdb, parquet, or pairs.")

    def partition_by(self, out_path: Path | str, partition_by: List[str] | str):
        partition_by = [partition_by] if isinstance(partition_by, str) else partition_by
        partition_by_str = ", ".join(partition_by)
        self.conn.execute(
f"""
COPY pairs TO '{out_path}' (
    FORMAT PARQUET, 
    PARTITION_BY ({partition_by_str}), 
    OVERWRITE_OR_IGNORE,
    KV_METADATA {{
        header: '{self.header.non_columns_text}'
    }}
)
"""
        )

    def iter_chroms(self, vector_multiple: int = None):
        """Yield Pandas dataframes in chunks in order of (chrom1, chrom2)
    
        Arguments:
            vector_multiple: The number of rows returned in each chunk is the vector size (2048 by default) * vector_multiple (1 by default).
        """
        combinations = self.conn.execute("SELECT DISTINCT chrom1, chrom2 FROM pairs ORDER BY chrom1, chrom2").pl()

        for chrom1, chrom2 in combinations.iter_rows():
            query = self.conn.execute(
                """
                SELECT * FROM pairs
                WHERE
                    chrom1 = $chrom1 
                    AND chrom2 = $chrom2
                ORDER BY chrom1, chrom2, pos1, pos2
                """,
                {"chrom1": chrom1, "chrom2": chrom2}
                            )
            if vector_multiple is None:
                yield query.df()
            else:
                while (chunk := query.fetch_df_chunk(vector_multiple)) is not None:
                    if chunk.empty:
                        break
                    else:
                        yield chunk
    

    @property
    def header(self) -> PairsHeader:
        header_text = self.conn.execute("SELECT * FROM metadata").fetchone()[1]
        return PairsHeader(header_text, from_text = True, columns = self.columns, replace_columns = True)

    @property
    def columns(self) -> List[str]:
        return [col[1] for col in self.conn.execute("PRAGMA table_info(pairs)").fetchall()]

    @property
    def conn(self) -> duckdb.DuckDBPyConnection:
        if self._conn is None:
            self._conn = duckdb.connect(self.path)
            if self.memory_limit:
                self._conn.set_memory_limit()
        return self._conn
    
    @property
    def current_database(self) -> str:
        return self.conn.execute("SELECT current_database()").fetchone()[0]
    
    def connect(self, out_path: Path | str = ":memory:") -> duckdb.DuckDBPyConnection:
        try:
            self._conn.close()
        except:
            pass
        self._conn = duckdb.connect(out_path)
        return self._conn

    @classmethod    
    def open(
            cls, 
            in_path: Path | str,
            in_format: str = "autodetect"
    ) -> "PairsSQL":
        """Open pairs, pairs, or parquet as PairsSQL object"""
        in_path = Path(in_path)
        if not in_format or in_format == "autodetect":
            in_format, _ = autodetect_extension(in_path)

        if in_format == "duckdb":
            result = PairsSQL().from_duckdb(in_path)
        elif in_format == "parquet":
            result = PairsSQL().from_parquet(in_path)
        elif in_format == "pairs":
            result = PairsSQL().from_pairs(in_path)
        else:
            raise Exception(f"Format {in_format} not recognized. Should be one of duckdb, parquet, or pairs.")


        return result