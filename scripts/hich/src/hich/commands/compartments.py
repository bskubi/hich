import click
from pathlib import Path

@click.command
@click.option("--chroms", multiple = True, default = None)
@click.option("--exclude-chroms", multiple = True, default = None)
@click.option("--keep-chroms-when", type = str, default = None)
@click.option("--n_eigs", type = int, default = 1)
@click.argument("reference")
@click.argument("matrix")
@click.argument("resolution", type = int)
def compartments(chroms, exclude_chroms, keep_chroms_when, n_eigs, reference, matrix, resolution):
    """
    """
    pass