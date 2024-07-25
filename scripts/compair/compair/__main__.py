import click
from smart_open import smart_open
from itertools import combinations
@click.command
@click.argument("f1")
@click.argument("f2")
def main(f1, f2):
    f1 = smart_open(f1, "rt").readlines()
    f2 = smart_open(f2, "rt").readlines()
    for i, lines in enumerate(zip(f1, f2)):
        l1, l2 = lines
        is_comment = l1.startswith("#") and l2.startswith("#")
        assert is_comment or l1 == l2, f"\nLine {i}\n{l1}{l2}"


if __name__ == "__main__":
    main()