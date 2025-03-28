import click
from hich.commands.pairs import pairs
from hich.commands.fasta import fasta
from hich.commands.report import report
from hich.commands.matrix import matrix

@click.group()
def hich():
    pass

hich.add_command(fasta)
hich.add_command(pairs)
hich.add_command(report)
hich.add_command(matrix)

if __name__ == "__main__":
    hich()