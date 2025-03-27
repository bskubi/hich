import click
from hich.commands.pairs import pairs
from hich.commands.fasta import fasta
from hich.commands.pipeline import pipeline
from hich.commands.matrix import matrix

@click.group()
def hich():
    pass

hich.add_command(fasta)
hich.add_command(pairs)
hich.add_command(pipeline)
hich.add_command(matrix)

if __name__ == "__main__":
    hich()