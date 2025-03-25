import click
from hich.commands.pairs import pairs

@click.group()
def hich():
    pass

hich.add_command(pairs)

if __name__ == "__main__":
    hich()