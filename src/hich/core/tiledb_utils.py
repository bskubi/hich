"""
Convenience methods for interacting with TileDB arrays.

Attributes:
    create (function): Create, update, or replace TileDb array.
    dim_names (function): Get names of TileDB array dims.
    attr_names (function): Get names of TileDB array attrs.
    write_pl_to_tiledb (function): Write a polars DataFrame to TileDB array.
"""
import tiledb
import shutil
import polars as pl
from pathlib import Path

def create(uri: str, schema: tiledb.ArraySchema, append: bool):
    """
    Remove any file/dir at URI, then create new SparseArray.

    It will be deleted if:
    1. 'append' is True and 'uri' is a TileDB array whose schema does not match.
    2. The target at 'uri' is not a TileDB array.
    3. 'append' is False.
    """
    target_exists = Path(uri).exists()

    if not target_exists:
        if not append:
            tiledb.SparseArray.create(uri, schema=schema)
        return

    # Determine if destruction is necessary
    requires_delete = False
    if not append:
        requires_delete = True
    else:
        try:
            requires_delete = tiledb.SparseArray(uri).schema != schema
        except:
            requires_delete = True

    # Execute destruction and recreation
    if requires_delete:        
        shutil.rmtree(uri, ignore_errors=True)
        tiledb.SparseArray.create(uri, schema=schema)

def dim_names(array: tiledb.Array) -> list[str]:
    """
    Get names of all TileDB Dims
    """
    return [array.domain.dim(i).name for i in range(array.domain.ndim)]

def attr_names(array: tiledb.Array) -> list[str]:
    """
    Get names of all tileDB Attrs
    """
    return array.attr_names

def write_pl_to_tiledb(
        df: pl.DataFrame,
        array: tiledb.SparseArray
    ):
    """
    Write polars DataFrame to TileDB (df must contain col names matching array names)
    """

    # Get matching column names from df
    _dim_names = dim_names(array)
    _attr_names = attr_names(array)

    # Reshape the DFs as expected by tiledb
    dims = tuple(df.select(*_dim_names).to_numpy().T)
    attrs = df.select(*_attr_names).to_dict() if _attr_names else None
    
    # Write to tiledb array
    array[dims] =  attrs