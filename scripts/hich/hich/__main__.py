from collections import defaultdict
from hich.cli import BooleanList, IntList, PathList, StrList
from hich.compartments import write_compartment_scores
from hich.fragtag import tag_restriction_fragments
from hich.hicrep_combos import hicrep_combos
from hich.visuals import view_hicrep
from pathlib import Path
from polars import DataFrame
import click
import hich.digest as _digest
import polars as pl
from polars import DataFrame
from hich.pairs import PairsClassifier, PairsFile
from hich.sample import SelectionSampler
from hich.stats import DiscreteDistribution, compute_pairs_stats_on_path_list, aggregate_classifier
from numbers import Number
from itertools import combinations

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

@hich.command
@click.option("--conjuncts",
    type = str,
    default = "record.chr1 record.chr2 record.pair_type stratum",
    show_default = True,
    help = "PairsSegment traits that define the category for each record (space-separated string list)")
@click.option("--cis_strata",
    type = str,
    default = "10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000",
    show_default = True,
    help = "PairsSegment cis distance strata boundaries (space-separated string list)")
@click.option("--orig_stats",
    type = str,
    default = "",
    show_default = True,
    help = ("Stats file containing original count distribution. Can be produced with hich stats. "
            "Computed from conjuncts and cis_strata if not supplied. Overrides default conjuncts and cis_strata if they are supplied."))
@click.option("--target_stats",
    type = str,
    default = "",
    show_default = True,
    help = "Stats file containing target count distribution.")
@click.option("--to_count",
    type = str,
    default = "",
    show_default = True,
    help = ("Float on [0.0, 1.0] for fraction of records to sample, or positive integer number of counts to sample. "
            "If a target stats file is supplied, further downsamples it to the given count."))
@click.argument("input_pairs_path", type = str)
@click.argument("output_pairs_path", type = str)
def downsample(conjuncts, cis_strata, orig_stats, target_stats, to_count, input_pairs_path, output_pairs_path):
    orig_classifier, orig_distribution = load_stats_and_classifier_from_file(orig_stats) if orig_stats else (None, None)
    target_classifier, target_distribution = load_stats_and_classifier_from_file(target_stats) if target_stats else (None, None)
    
    if orig_classifier and target_classifier:
        assert orig_classifier.conjuncts == target_classifier.conjuncts, f"Original and target conjuncts do not match for {orig_stats} and {target_stats}."
        conjuncts = orig_classifier.conjuncts
        if "stratum" in orig_classifier and "stratum" in target_classifier:
            cis_strata = list(set(orig_classifier["stratum"] + target_classifier["stratum"]))        
    elif orig_classifier:
        conjuncts = orig_classifier.conjuncts
    elif target_classifier:
        conjuncts = target_classifier.conjuncts
    
    classifier = PairsClassifier(conjuncts, cis_strata)
    
    if to_count.isdigit():
        to_count = int(to_count)
    else:
        try:
            to_count = float(to_count)
        except ValueError:
            to_count = None
    
    if not orig_distribution:
        _, orig_distribution = compute_pairs_stats_on_path((classifier, pairs_path))
    if not target_distribution:
        assert to_count is not None, "No target distribution or count supplied for downsampling."
        target_distribution = orig_distribution.to_count(to_count)
    
    sampler = SelectionSampler(orig_distribution, target_distribution)
    input_pairs_file = PairsFile(input_pairs_path)
    output_pairs_file = PairsFile(output_pairs_path, mode = "w", header = input_pairs_file.header)

    for record in input_pairs_file:
        outcome = classifier.classify(record)
        if sampler.sample(outcome):
            output_pairs_file.write(record)

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
@click.option("--format", "--fmt", "fmt",
    type = click.Choice(['fasta', 'fastq', 'seqio', 'sam', 'bam', 'sambam', 'alignment', 'pairs']), required = True,
    help = "Alignment data format")
@click.option("--f1", "-f", "--file", "--file1", type = str, required = True,
    help = "Path to first (or only) input sequencing data file")
@click.option("--f2", "--file2", "f2", default = None, type = str, show_default = True,
    help = "Path to second input sequencing data file")
@click.option("--out-dir", "--dir", "--output-dir", type = str, default = "",
    help = "Output directory")
@click.option("--annot-file", "--annot", "-a", "--annotations", default = None, type = str,
    help = ("Path to annotation file, a columnar text file mapping data file "
            "information such as a cell barcode to new information such as "
            "an experimental condition or cell ID. Annotation files with headers "
            "convert to a dict with format {col1_row: {col2:col2_row, col3:col3_row...}}."))
@click.option("--annot-has-header", "-h", type = bool, default = False, show_default = True,
    help = "Whether or not annotation file has a header row.")
@click.option("--annot-separator", "-s", type = str, default = "\t", show_default = repr('\t'),
    help = "The column separator character")
@click.option("--head", "--n_records", type = int, default = None, show_default = True,
    help = "Take only the first n records from the data files. Takes all records if not provided.")
@click.option("--key-code", "--kc", type = str, default = None,
    help = "Python code to extract an annotation row key from the current record.")
@click.option("--record-code", "--rc", type = str, default = None,
    help = "Python code to modify the record")
@click.option("--output-code", "--fc", "--filename", type = str, default = None,
    help = "Python code to select the record output.")
def organize(fmt,
             f1,
             f2,
             out_dir,
             annot_file,
             annot_has_header,
             annot_separator,
             head,
             key_code,
             record_code,
             output_code):
    """Reannotate and split sequencing data files


    """
    raise NotImplementedError("Hich organize is not implemented yet")

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

@hich.command
@click.option("--conjuncts",
    type = str,
    default = "record.chr1 record.chr2 record.pair_type stratum",
    show_default = True,
    help = "PairsSegment traits that define the category for each record (space-separated string list)")
@click.option("--cis_strata",
    type = str,
    default = "10 20 50 100 200 500 1000 2000 5000 10000 20000 50000 100000 200000 500000 1000000 2000000 5000000",
    show_default = True,
    help = "PairsSegment cis distance strata boundaries (space-separated string list)")
@click.option("--stats_directory",
    type = str,
    default = "",
    show_default = True,
    help = "Stats tab-separated file saved at [stats_directory]/[pairs_path].stats.tsv")
@click.argument("pairs_paths",
                type = str,
                nargs = -1)
def stats(conjuncts, cis_strata, stats_directory, pairs_paths):
    assert pairs_paths, "No pairs paths supplied."
    conjuncts = conjuncts.split()
    split_strata = lambda: [int(stratum) for stratum in cis_strata.split()]
    cis_strata = split_strata() if cis_strata else None
    classifier = PairsClassifier(conjuncts, cis_strata)

    result = compute_pairs_stats_on_path_list(classifier, pairs_paths)
    
    for pairs_path, distribution in result:
        df = classifier.to_polars(distribution)
        output_path = str(Path(stats_directory) / f"{pairs_path}.stats.tsv")
        print(output_path)
        df.write_csv(output_path, separator="\t", include_header=True)


@hich.command
@click.option("--to_group_mean", is_flag = True, default = False)
@click.option("--to_group_min", is_flag = True, default = False)
@click.option("--to_count", type = str, default = None)
@click.option("--prefix", type = str, default = None)
@click.argument("stats_paths", type = str, nargs = -1)
def stats_aggregate(to_group_mean, to_group_min, to_count, prefix, stats_paths):
    """Aggregate hich stats files for .pairs

    """
    classifier, distributions = aggregate_classifier(stats_paths)
    targets = [d for d in distributions]
    build_prefix = ""
    if to_group_mean:
        group_mean = DiscreteDistribution.mean_mass(distributions)
        targets = [d.downsample_to_probabilities(group_mean) for d in distributions]
        if prefix is None:
            build_prefix += "to_group_mean"
    if to_group_min:
        total = min(targets)
        targets = [d.to_count(total) for d in targets]
        if prefix is None:
            build_prefix += "to_group_min"
    if to_count:
        if to_count.isdigit():
            to_count = int(to_count)
        else:
            to_count = float(to_count)
        targets = [d.to_count(to_count) for d in targets]
        if prefix is None:
            build_prefix += f"to_{to_count}"
    if prefix is None:
        prefix = build_prefix + "_"
    for d, stats_path in zip(targets, stats_paths):
        df = classifier.to_polars(d)
        path = str(Path(stats_path).parent / (prefix + Path(stats_path).name))
        print(path)
        df.write_csv(path, separator = "\t", include_header = True)
    



    




@hich.group 
def view(): pass

@view.command(name='hicrep')
@click.option("--host", default = "127.0.0.1", show_default = True)
@click.option("--port", default = 8050, show_default = True)
def hicrep_comparisons(host, port):
    view_hicrep.run_dashboard(host, port)

if __name__ == "__main__":
    hich()
