import click
import hictkpy
from hich.matrix.hicrep import HiCRepCaller
import itertools
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import smart_open_with_pbgzip
from smart_open import smart_open

@click.group
def matrix():
    pass

@matrix.command
@click.option("--resolution", "-r", "resolutions", type=int, multiple=True, 
              help = "Contact matrix resolutions")
@click.option("--h", "-h", multiple=True, type=int, default = [1], 
              help = "Values of hicrep smoothing (h) parameters")
@click.option("--d-bp-max", "--dBPMax", "-m", multiple=True, type = int, default = [-1], 
              help = "Values of hicrep dBPMax parameter")
@click.option("--b-downsample", "--bDownSample", "-d", multiple=True, type=bool, default = [False], 
              help = "Values of hicrep bDownSample parameter")
@click.option("--normalization", "-n", "normalizations", multiple=True, type=str, default = None,
              help = "Contact matrix normalizations. If unspecified, uses all available; set --raw to include un-normalized with all the other normalizations. Use 'NONE' to call only on raw contact matrix.")
@click.option("--raw", is_flag = True, default = False, 
              help = "Call HiCRep without normalization (normalization = 'NONE') in addition to any other normalizations specified.")
@click.option("--one-norm", is_flag = True, default = False, help = "If multiple normalizations are present, only use one of them.")
@click.option("--chrom", "-c", "chroms", multiple=True, type=str, default=None, help = "The set of chromosomes to use. If unspecified, SCC scores are called on the subset of chromosomes shared by all contact matrices.")
@click.option("--skip-chrom", "-s", "skip_chroms", multiple=True, type=str, default=None, help = "Exclude these chromosomes, even if specified by --chrom")
@click.option("--region", "-g", multiple = True, type = str, default = None, help = "UCSC genome region to compare, i.e. chr1:10000000-11000000")
@click.option("--bed", "-b", multiple = True, type = str, default = None, help = "BED file defining a partition over the genome with HiCRep called on each region")
@click.option("--partition", "-p", multiple = True, type = int, default = None, help = "Partition selected chromosomes into uniform regions of this size")
@click.option("--n_proc", type=int, default=None, help = "Number of parallel processes to use to call SCC scores. By default, uses all available processors.")
@click.argument("out_path", type=str)
@click.argument("in_paths", type=str, nargs = -1)
def hicrep(resolutions, h, d_bp_max, b_downsample, normalizations, raw, one_norm, chroms, skip_chroms, region, bed, partition, n_proc, out_path, in_paths):
    """Call HiCRep SCC scores on combinations of contact matrices and parameterizations
    
    \b
    OUT_PATH: Output filename for results. The caller parameters are stored along with the resulting SCC score in the order returned.
    IN_PATHS: Two or more .hic or .cool/.mcool files (does not work on .scool files)

    \b
    Examples:

    \b
    Call SCC scores on all chromosomes shared by all contact matrices at 10kb and 100kb resolution using h values of 0, 1 and 2, scoring results in scc.txt:
        hich matrix hicrep -r 10000 -r 100000 -h 0 -h 1 -h 2 scc.txt mx1.mcool mx2.mcool mx3.mcool
    """
    # Get common chromosomes and normalizations for all files
    in_paths = set([Path(in_path).absolute() for in_path in in_paths])
    
    # Validate that we have at least two files and that all files exist
    assert all([it.exists() for it in in_paths]), f"Some input contact matrices do not exist on disk: {[it for it in in_paths if not it.exists()]}" 
    assert len(in_paths) >= 2, f"At least two contact matrices are required, but only {len(in_paths)} was suppled. IN_PATHS was {in_paths}."
    in_paths = [str(it) for it in in_paths]

    # Update chromosomes and normalizations
    chroms = set(chroms)
    normalizations = set(normalizations)

    if not chroms or not normalizations:
        # Get chromosomes and normalizations common to all datasets if not specified by the user
        update_chroms = None
        update_normalizations = None
        for path, resolution in itertools.product(in_paths, resolutions):
            file = hictkpy.File(path, resolution)
            if not chroms:
                new_chroms = set(file.chromosomes().keys())
                update_chroms = new_chroms if update_chroms is None else update_chroms.intersection(new_chroms)
            if not normalizations:
                new_normalizations = set(file.avail_normalizations())
                update_normalizations = new_normalizations if new_normalizations is None else new_normalizations.intersection(new_normalizations)
        chroms = update_chroms if not chroms else chroms
        normalizations = update_normalizations if not normalizations else normalizations

    # Add directive to use un-normalized matrix if requested
    if raw:
        normalizations.add("NONE")

    # Select a single normalization if requested
    if one_norm:
        chroms = set(list(chroms)[0])

    # Remove skipped chromosomes
    for chrom in skip_chroms:
        if chrom in chroms:
            chroms.discard(chrom)

    # Make sure we still have at least one chromosome and normalization
    assert chroms, f"No shared chromosomes amongst input matrices {in_paths}"
    assert normalizations, f"No shared normalizations amongst input matrices {in_paths}"

    # Create all individual combinations of parameters
    paths = itertools.combinations_with_replacement(in_paths, r=2)
    params = itertools.product(paths, chroms, normalizations, resolutions, h, d_bp_max, b_downsample)
    callers = [HiCRepCaller(*caller[0], *caller[1:]) for caller in params]

    # Call HiCRep SCC scores and write as TSV output
    f = smart_open(out_path, "wt")

    with ProcessPoolExecutor(max_workers=n_proc) as ppe:
        header_written = False
        for result in ppe.map(HiCRepCaller.call_scc, callers):
            if not header_written:
                header = "\t".join(result.as_dict().keys()) + "\n"
                f.write(header)
            row = "\t".join([str(val) for val in result.as_dict().values()]) + "\n"
            f.write(row)
    f.close()