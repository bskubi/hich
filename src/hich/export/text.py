import tiledb
import click
import sys
import duckdb
from dataclasses import dataclass, field
import pandas as pd
from typing import Iterator

@dataclass
class ExportTextEngine:
    # TileDB array to export from
    array: tiledb.SparseArray

    # Filter and process
    query_order: str
    cond: str | None
    sql: str | None

    # Output formatting
    sep: str
    header: bool
    dims_attrs: list[str] = field(default_factory=list)
    index_label: str | None = None
    index: bool = False

    def use_format(self, format: str):
        match format:
            case "4dn-pairs":
                self.index=True
                self.index_label = "readID"
                self.dims_attrs = ["chrom1", "pos1", "chrom2", "pos2", "strand1", "strand2"] + self.dims_attrs

                if "readID" in self.all_names:
                    self.dims_attrs = ["readID"] + self.dims_attrs
                    self.index = False
                    self.index_label = None


            case "tsv":
                self.dims_attrs = list(self.dims_attrs)
                self.index_label=None
                self.index=False
            case _:
                raise NotImplementedError(f"Format {format} not implemented.")

    @property
    def all_dim_names(self) -> list[str]:
        return [self.array.dim(i) for i in range(self.array.ndim)]

    @property
    def all_attr_names(self) -> list[str]:
        return self.array.attr_names
    
    @property
    def all_names(self) -> list[str]:
        return self.all_dim_names + self.all_attr_names
    
    @property
    def dims(self) -> list[str]:
        return [dim for dim in self.all_dim_names if dim in self.dims_attrs]
    
    @property
    def attrs(self) -> list[str]:
        return [attr for attr in self.all_attr_names if attr in self.dims_attrs]
    
    @property
    def sep_unicode(self) -> str:
        return self.sep.encode().decode('unicode_escape')

    @property
    def query_order_tiledb(self) -> str:
        match self.query_order:
            case "global": return  "G"
            case "col-major": return "F"
            case "row-major": return "C"
            case "unordered": return "U"

    @property
    def query(self) -> tiledb.Query:
        return self.array.query(attrs=self.attrs, cond=self.cond or None, dims=self.dims, order=self.query_order_tiledb, return_incomplete=True)
    
    @property
    def chunks(self) -> Iterator[tuple[int, pd.DataFrame]]:
        for i, df_chunk in enumerate(self.query.df[:]):
            # Results are chunked for memory safety. Extract required dims and attrs
            if self.dims_attrs:
                df_filtered = df_chunk[self.dims_attrs]
            else:
                df_filtered = df_chunk
            yield i, df_filtered

    def export_text_sql(self):
        with duckdb.connect() as con:
            for i, df in self.chunks:
                if i == 0:
                    con.execute("CREATE TABLE query AS SELECT * FROM df")
                else:
                    con.execute("INSERT INTO query SELECT * FROM df")
            con.execute(f"COPY ({self.sql}) TO '/dev/stdout' (FORMAT CSV, DELIMITER '{self.sep_unicode}', HEADER {str(self.header).lower()})")

    def export_text_direct(self):
        for i, df in self.chunks:
            df.to_csv(sys.stdout, sep=self.sep_unicode, header=self.header, index_label=self.index_label, index=self.index)

    def export_text(self):
        self.export_text_sql() if self.sql else self.export_text_direct()

@click.command()
@click.option(
    "--format", 
    type=click.Choice(["tsv", "4dn-pairs"]), 
    default="tsv",
    help="Format of output matrix."
)
@click.option("--query-order", type=click.Choice(["row-major", "col-major", "global", "unordered"]), default="row-major", help="Sort order of rows extracted from TileDB")
@click.option("--cond", type=str, default=None, help="Selection condition.")
@click.option("--sql", type=str, default=None, help="Store selected rows in intermediate DuckDB table 'query', then export results of SQL command.")
@click.option("--sep", default=r"\t", help="Column separator")

@click.option("--header", is_flag=True, default=False, help="Output with header")
@click.argument("path")
@click.argument("dims_attrs", nargs=-1)
def export_text(format, query_order, cond, sql, header, sep, path, dims_attrs):
    """
    Export rows from TileDB with tab-delimited columns.
    """
    engine = ExportTextEngine(
        array = tiledb.SparseArray(path),
        query_order=query_order,
        cond=cond,
        sql=sql,
        sep=sep,
        header=header,
        dims_attrs=list(dims_attrs)
    )
    engine.use_format(format)
    engine.export_text()


if __name__ == "__main__":
    export_text()