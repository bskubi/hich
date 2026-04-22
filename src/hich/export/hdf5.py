import click
import tiledb
import h5py
from .query_engine import DuckDBExportEngine, ExportEngine, HDF5WriteStrategy, CSVWriteStrategy, ParquetWriteStrategy, SciPySparseWriteStrategy

@click.command()
@click.option("--query-name", "--dim", "--attr", "--name", "-n", "query_names", multiple=True)
@click.option("--query-cond", type=str, default=None)
@click.option("--query-order", type=click.Choice(["global", "col-major", "row-major", "unordered"]), default="col-major")
@click.option("--sql", multiple=True)
@click.option("--metadata", type=(str, str, str), multiple=True)
@click.option("--format", type=click.Choice(["hdf5", "csv", "parquet", "csr", "amethyst.h5", "bismark.cov"]))
@click.option("--dtype", type=(str, str), multiple=True)
@click.option("--header", is_flag=True, default=False)
@click.option("--sep", type=str, default="\t")
@click.option("--ext", type=str, default=None)
@click.option("--partition-by", type=str, default="")
@click.argument("tiledb_path")
@click.argument("output_path")
def export_partition(query_names, query_cond, query_order, sql, metadata, format, dtype, header, sep, ext, partition_by, tiledb_path, output_path):
    """
    Export from TileDB to a partition on disk or within an HDF5 file.
    """
    partition_by = partition_by.strip("/").split("/")
    sql = list(sql)
    metadata = list(metadata)
    custom_dtypes = dict(dtype) if dtype else None
    query_names = list(query_names) if query_names else None

    match format:
        case "amethyst.h5":
            # Amethyst single bp resolution datasets
            query_names = ["context", "barcode", "chrom", "pos", "mC"]
            custom_dtypes = {"chrom": "S20", "pos": "int64"}
            partition_by = ["context", "barcode", "1"]

            sql=[
                "SELECT context, barcode, chrom, pos, SUM(mC=0) AS t, SUM(mC=1) AS c "
                "FROM QUERY_RESULT "
                "GROUP BY context, barcode, chrom, pos " 
                "ORDER BY context, barcode, chrom, pos"
            ]
            metadata=[("/metadata/version", "amethyst2.0.0", "S20")]
            strategy = HDF5WriteStrategy(h5py.File(output_path, "a"), custom_dtypes=custom_dtypes)
        case "bismark.cov":
            # Output from MethSCan prepare
            # CSR matrix (pos x cells)
            # shape: (chrom_size+1, n_cells)
            # dtype: np.int8
            # values: 1: methylated, -1 unmethylated
            query_names = ["barcode", "chrom", "pos", "mC"]
            custom_dtypes = {"chrom": "S20"}
            partition_by = ["barcode"]
            sql = [
                "SELECT barcode, chrom, pos AS 'start', pos+1 AS 'end', SUM(mC=1)/(COUNT(*))*100 AS pct, SUM(mC=0) AS t, SUM(mC=1) AS c "
                "FROM QUERY_RESULT "
                "GROUP BY barcode, chrom, pos " 
                "ORDER BY barcode, chrom, pos"
            ]
            strategy = CSVWriteStrategy(output_path, header=False, sep="\t", ext=".cov")
        case "hdf5":
            strategy = HDF5WriteStrategy(h5py.File(output_path, "a"), custom_dtypes=custom_dtypes, ext=ext or ".h5")
        case "csv":
            strategy = CSVWriteStrategy(output_path, header=header, sep=sep, ext=ext)
        case "parquet":
            strategy = ParquetWriteStrategy(output_path, ext=ext or ".parquet")
        case "csr":
            strategy = SciPySparseWriteStrategy(output_path, ext=ext or ".npz")
        case _: raise NotImplementedError(f"Format {format} not implemented.")

    tiledb_array = tiledb.SparseArray(tiledb_path)

    engine_kwargs = {"tiledb_array": tiledb_array, 
                     "query_names": query_names, 
                     "query_cond": query_cond,
                     "query_order": query_order}
    if sql:
        engine = DuckDBExportEngine(sql=sql, **engine_kwargs)
    else:
        engine = ExportEngine(**engine_kwargs)
    
    engine.write_partition(strategy, partition_by)
    
    # for _metadata in metadata:
    #     file.create_dataset(_metadata[0], data=_metadata[1], dtype=_metadata[2])