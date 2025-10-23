from smart_open import smart_open
import smart_open_with_pbgzip
import click

@click.command
@click.argument("fq1")
@click.argument("fq2")
@click.argument("fq-out")
def interleave(fq1, fq2, fq_out):
    try:
        lines1 = smart_open(fq1, "rt").readlines()
        lines2 = smart_open(fq2, "rt").readlines()

        count = len(lines1) + len(lines2)
        output = smart_open(fq_out, "w")
        l1 = 0
        l2 = 0
        while l1 < len(lines1) or l2 < len(lines2):
            if l1 <= l2:
                lines = lines1[l1:l1+4]
                l1 += 4
            else:
                lines = lines2[l2:l2+4]
                l2 += 4
            for line in lines:
                output.write(line)
    except Exception as e:
        raise Exception(f"With arguments {fq1}, {fq2}, {fq_out}") from e
if __name__ == "__main__":
    interleave()


