from hich.parse.pairs_file import PairsFile
from hich.coverage.pairs_space import PairsSpace
from hich.coverage.trans_cis_thresholds import TransCisThresholds
from hich.coverage.pairs_histogram import PairsHistogram
import click
from dataclasses import dataclass, field
from collections import defaultdict

def test():
    pairs_file = PairsFile("data/240802_Arima_reseq_dedup.pairs.gz")
    cis_thresholds = [10, 20, 50, 100]
    space = TransCisThresholds(ignore_code = "not pair.ur()", cis_thresholds = cis_thresholds)

    dist = PairsHistogram(space)

    for i, pair in enumerate(pairs_file):
        dist.count_pair(pair)
        if i >= 1000:
            break

    print(dist)
    print(dist.to_probdist())
    prob = [.7, .15, .05, .1, 0]
    print(dist.downsample_to_probdist(prob))
    print(dist.downsample_to_probdist(prob).to_count(50))

from json import dumps, loads
from pathlib import Path
import polars as pl
import re
default_cis_thresholds = [10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000]

def extract_integers(input_string):
    # Use regex to find all sequences of digits
    return [int(num) for num in re.findall(r'\d+', input_string)]

class DownsampleOption(click.ParamType):
    name = "downsample_option"

    def convert(self, value, param, ctx):
        # Check for allowed strings first
        if value in ["min_in_group", "min_all_groups"]:
            return value

        # Check for a number between 0 and 1 (inclusive)
        if re.fullmatch(r'0(\.\d+)?|1(\.0+)?', value):
            return float(value)
        
        # Check for an integer 1 or greater
        if re.fullmatch(r'[1-9]\d*', value):
            return int(value)

        msg = "\n".join([
            f"{value} is not a valid option for --downsample. Must be chosen from:",
            "\t- [0.0-1.0]: Exact fraction of read count",
            "\t- 1+, integer: Exact number of final reads",
            "\t- min_in_group: Downsample to size of smallest sample in the read group",
            "\t- min_all_groups: Downsample to size of smallest sample across all read groups"
        ])
        self.fail(msg, param, ctx)

def load_groups(groups: str, sep: str = "\t", list_sep: str = ",") -> dict[str, str]:
    # This should first consider the groups as a list of pairs
    # if not, then a filename for a csv
    # if not, a JSON string
    if Path(groups).exists():
        pass
        #groups = pl.read_csv(groups, separator = sep, has_header = False).to_dict()
    else:
        try:
            groups = loads(groups)
        except:
            start = f"--groups is '{groups}', but must be either of:"
            csv = f"{repr(sep)}-delimited headerless two-column file with files in column 1 and groups in column 2"
            json = "JSON dict with filenames as keys and group ids as values."
            msg = start + f"\n\t- {csv}" + f"\n\t- {json}"
            print(msg)
            exit()
    return groups

def validate_filenames(filenames: list[str]):
    for file in filenames:
        assert Path(file).exists(), f"File {file} not found."

def read_counts_histogram(data: tuple[PairsFile, PairsHistogram, int]) -> PairsHistogram:
    pairs_file, histogram, lines = data
    for i, pair in enumerate(pairs_file):
        histogram.count_pair(pair)
        if lines is not None and i >= lines:
            break
    return histogram

@dataclass
class PairsFileData:
    pairs_file: PairsFile
    histogram: PairsHistogram

@dataclass
class PairsFileDataGrouping:
    groups: defaultdict(list) = field(default_factory = lambda: defaultdict(list))
    ignore: list[str] = field(default_factory = list)

    def to_centers(self):
        for group, data in self.groups.items():
            use_histograms = [d.histogram for d in data if d.pairs_file.filepath_or_object not in ignore]
            central_histogram = PairsHistogram.center(use_histograms)
            for i, datum in enumerate(data):
                self.groups[group][i].histogram = central_histogram
    
    def 


@click.command
@click.option("--groups", "-g", "--sample_groups", type = str, default = {})
@click.option("--downsample", type=DownsampleOption(), default = "min_all_groups")
@click.option("--outlier", multiple = True)
@click.option("--ignore-contig", "--ic", multiple = True)
@click.option("--ignore-cis", type = bool, default = False, show_default = True)
@click.option("--ignore-trans", type = bool, default = False, show_default = True)
@click.option("--cis-thresholds", "--thresholds", "--cis-strata", "--strata", type = str, default = default_cis_thresholds, show_default = True)
@click.option("--sep", default = "\t")
@click.option("--list-sep", default = ",")
def coverage(groups, downsample, outlier, ignore_contig, ignore_cis, ignore_trans, cis_thresholds, sep, list_sep):
    groups = load_groups(groups, sep, list_sep)
    if isinstance(groups, str):
        groups = {groups: 0}
    cis_thresholds = extract_integers(cis_thresholds)
    space = TransCisThresholds(ignore_code = "not pair.ur()", cis_thresholds = cis_thresholds)
    
    data_grouping = PairsFileDataGrouping()
    print(groups)
    for file, group in groups.items():
        pairs_file = PairsFile(file)
        histogram = PairsHistogram(space)
        histogram = read_counts_histogram((pairs_file, histogram, 10000))
        data_grouping.groups[group] = PairsFileData(pairs_file, histogram)
    
    #histogram.write_csv("hist.tsv", separator="\t")
    #validate_filenames(groups.keys())

if __name__ == "__main__":
    coverage()

"""
0. On input
    Associate samples with downsampling groups
        Default: all in their own target profile group

    Specify a minimum number of reads per contig.
        - Above: follow normal methods
        - Below: either retain all or drop all

1. To form target probdist PairHistograms, convert count PairHistograms to
probabilities and average them, leaving profiles unchanged if they are alone
in the sample group.
    - Can also provide an explicit target for each sample group
    - Specify some sample files and/or contigs to ignore

2. The user can specify:
    2a) (optional) Downsample to target probdist for sample group
    2b) (optional) Downsample to any of:
        - Specific size
        - Specific fraction
        - Minimum sample size within sample group
        - Minimum sample size across sample groups

3. This yields new PairsHistograms for each sample. We then use MultiSampler
to selection sample the pairs files.

Possibly reading from the PairsHistograms should be parallelized
"""