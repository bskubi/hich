import click
from hich.cli import BooleanList, IntList, PathList, StrList
import hich.digest as _digest
from hich.fragtag import tag_restriction_fragments
from hich.hicrep_combos import hicrep_combos
import hich.coverage as _coverage
from hich.compartments import write_compartment_scores
import polars as pl
from pathlib import Path
from hich.organize import organize as _organize

@click.group
def hich():
    pass

@hich.command
@click.option("--chroms", type = StrList, default = None)
@click.option("--exclude-chroms", type = StrList, default = None)
@click.option("--keep-chroms-when", type = str, default = None)
@click.option("--n_eigs", type = int, default = 1)
@click.argument("reference")
@click.argument("matrix")
@click.argument("resolution", type = int)
def compartments(chroms, exclude_chroms, keep_chroms_when, n_eigs, reference, matrix, resolution):
    matrix = Path(matrix)
    reference = Path(reference)
    final_suffix = matrix.suffixes[-1]
    prefix = matrix.name.rstrip(final_suffix)

    write_compartment_scores(prefix, matrix, reference, resolution, chroms, exclude_chroms, keep_chroms_when, n_eigs)

@hich.command()
@click.option("--strata", type = IntList, default = "0")
@click.option("--fraction", type = float, default = 1.0)
@click.argument("filenames", nargs = -1)
def coverage(strata, fraction, filenames):
    import warnings
    warnings.warn("Coverage not fully implemented yet.", UserWarning)

    strata = [0] if not strata else strata
    stratum_counts = _coverage.strata(filenames[0], strata)
    stratum_counts = stratum_counts.rename({filenames[0]:"N"})
    
    stratum_counts = stratum_counts.with_columns((pl.col("N")*fraction).cast(pl.Int64).alias("n"))
    _coverage.selection_sample(filenames[0], "output.pairs", strata, stratum_counts, ["chrom1", "chrom2", "stratum"], "N", "n")

@hich.command()
@click.option("--output", default = None, show_default = True, help = "Output file. Compression autodetected by file extension. If None, prints to stdout.")
@click.option("--startshift", default = 0, show_default = True, help = "Fixed distance to shift start of each fragment")
@click.option("--endshift", default = 0, show_default = True, help = "Fixed distance to shift end of each fragment")
@click.option("--cutshift", default = 1, show_default = True, help = "Fixed distance to shift cutsites")
@click.argument("reference")
@click.argument("digest", nargs = -1)
def digest(output, startshift, endshift, cutshift, reference, digest):
    """
    In silico digestion of a FASTA format reference genome into a
    BED format fragment index.

    Allows more than 800 restriction enzymes (all that are supported by
    Biopython's Restriction module, which draws on REBASE).

    Digest can also be specified as a kit name. Currently supported kits:

        "Arima Genome-Wide HiC+" or "Arima" -> DpnII, HinfI
    
    Multiple kits and a mix of kits and enzymes can be added. Duplicate
    kits, enzyme names, and restriction fragments are dropped.

    Can read compressed inputs using decompression autodetected and supported
    by the Python smart_open library. Can compress output using compression
    formats autodetected by the Polars write_csv function.

    Format is:
    chrom start cut_1
    chrom cut_1 cut_2

    The startshift param is added to all values in column 1.
    The endshift param is added to all values in column 2.
    """
    # We aim to support specification of digests by kit name
    # (potentially versioned), so this converts the kit names to the enzymes
    # used in that kit.
    _digest.digest(output, startshift, endshift, cutshift, reference, digest)

@hich.command
@click.option("--batch_size", default = 1000000)
@click.argument("fragfile")
@click.argument("out_pairs")
@click.argument("in_pairs")
def fragtag(batch_size, fragfile, out_pairs, in_pairs):
    tag_restriction_fragments(fragfile, in_pairs, out_pairs, batch_size)

@hich.command
@click.option("--format", "--fmt", "fmt", type = str, required = True)
@click.option("--f1", "-f", "--file", "--file1", type = str, required = True)
@click.option("--f2", "--file2", "f2", default = None, type = str)
@click.option("--out-dir", "--dir", "--output-dir", type = str, default = "")
@click.option("--annot-file", "--annot", "-a", "--annotations", default = None, type = str)
@click.option("--annot-has-header", "-h", type = bool, default = False, show_default = True)
@click.option("--annot-separator", "-s", type = str, default = "\t", show_default = True)
@click.option("--head", type = int, default = None, show_default = True)
@click.option("--record-code", "--rc", type = str, default = None)
@click.option("--filename-code", "--fc", "--filename", type = str, default = None)
def organize(fmt, f1, f2, out_dir, annot_file, annot_has_header, annot_separator, head, record_code, filename_code):
    _organize(fmt, f1, f2, out_dir, annot_file, annot_has_header, annot_separator, head, record_code, filename_code)

@hich.command
@click.option("--resolutions", type = IntList, default = 10000)
@click.option("--chroms", "--include_chroms", type = StrList, default = None)
@click.option("--exclude", "--exclude_chroms", type = StrList, default = None)
@click.option("--chromFilter", type=str, default = "chrom if size > 5000000 else None")
@click.option("--h", type = IntList, default = "1")
@click.option("--d_bp_max", type = IntList, default = "-1")
@click.option("--b_downsample", type = BooleanList, default = False)
@click.option("--nproc", type=int, default=None)
@click.option("--output", type=str, default = None)
@click.argument("paths", type=str, nargs = -1)
def hicrep(resolutions, chroms, exclude, chromFilter, h, d_bp_max, b_downsample, nproc, output, paths):
    result = hicrep_combos(resolutions, chroms, exclude, chromFilter, h, d_bp_max, b_downsample, nproc, output, paths)
    if result is not None:
        click.echo(result)
    


if __name__ == "__main__":
    hich()
