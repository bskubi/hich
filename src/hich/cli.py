import click
from .extract.extract_dna_methylation_engine import extract_dna_methylation
from .extract.extract_3d_pairs_engine import extract_3d_pairs
from .export.text import export_text
from .export.hdf5 import export_partition

@click.group()
def hich():
    pass


@hich.group()
def extract():
    pass

extract.add_command(extract_dna_methylation, name="dna-methylation")
extract.add_command(extract_3d_pairs, name="3d-pairs")

@hich.group()
def export():
    pass

export.add_command(export_text, name="text")
export.add_command(export_partition, name="partition")