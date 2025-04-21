import click
from hich.commands.pairs import pairs
from hich.commands.fasta import fasta
from hich.commands.workflow import workflow
from hich.commands.matrix import matrix

@click.group()
def hich():
    """CLI tools to process Hi-C data
    
    Version: 0.2.7
    """
    pass

hich.add_command(fasta)
hich.add_command(pairs)
hich.add_command(matrix)
hich.add_command(workflow)

if __name__ == "__main__":
    hich()
