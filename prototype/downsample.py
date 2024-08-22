from hich.parse.pairs_file import PairsFile
from hich.coverage.pairs_space import PairsSpace
from hich.coverage.trans_cis_thresholds import TransCisThresholds
from hich.coverage.pairs_histogram import PairsHistogram
from hich.sample.sampler import Sampler

import click
from dataclasses import dataclass, field
from collections import defaultdict
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

def read_counts_histogram(data: tuple[PairsFile, PairsHistogram, int]) -> PairsHistogram:
    pairs_file, histogram, lines = data
    for i, pair in enumerate(pairs_file):
        histogram.count_pair(pair)
        if lines is not None and i >= lines:
            break
    return histogram


@dataclass
class PairsData:
    filename: str
    outlier: bool
    pairs_file: PairsFile = None
    counts: PairsHistogram = None
    target_dist: PairsHistogram = None
    target_count: int = None
    samplers: defaultdict(Sampler) = None

    def compute_counts(self):
        for i, pair in enumerate(self.pairs_file):
            self.counts.count_pair(pair)
            if i > 1000:
                break
    
    def create_samplers(self):
        self.samplers = defaultdict(Sampler)
        for event in self.target_dist.events():
            orig_count = self.counts.distribution[event]
            target_count = self.target_dist.distribution[event]
            self.samplers[event] = Sampler(orig_count, target_count)
    
    def write_samples(self, write_file):
        read_from = PairsFile(self.filename)
        write_to = PairsFile(write_file, header = read_from.header, mode = "w")
        for pair in read_from:
            event = str(self.counts.space.event(pair))
            if self.samplers[event].sample(): write_to.write(pair)
            if all([sampler.finished() for sampler in self.samplers.values()]):
                break

    def __str__(self):
        s = []
        s.append(self.filename)
        s.append(f"\t- outlier: {self.outlier}")
        s.append(f"\t- counts: " + str(self.counts))
        s.append(f"\t- target_dist: " + str(self.target_dist))
        s.append(f"\t- target_count: " + str(self.target_count))
        return "\n".join(s)

@dataclass
class PairsComparison:
    space: PairsSpace
    groups: defaultdict[list[PairsData]]

    def __init__(self, space: PairsSpace, groups: dict[str: object], target_count: str = None):
        self.space = space
        self.groups = defaultdict(list[PairsData])
        for file, group in groups.items():
            data = PairsData(file, False, PairsFile(file), PairsHistogram(self.space), target_count = target_count)
            self.groups[str(group)].append(data)

    def compute_counts(self):
        for group, data in self.groups.items():
            for datum in data:
                datum.compute_counts()

    def compute_central_distributions(self):
        for group, data in self.groups.items():
            non_outliers = [datum.counts for datum in data if not datum.outlier]
            center = PairsHistogram.center(non_outliers)
            for i, datum in enumerate(data):
                self.groups[group][i].target_dist = center
    
    def downsample_to_target_probdist(self):
        for group, data in self.groups.items():
            for datum in data:
                datum.target_dist = datum.counts.downsample_to_probdist(datum.target_dist)
    
    def compute_target_count(self):
        min_in_groups = {}
        min_all_groups = 0
        positive_int = lambda v: isinstance(v, int) and v >= 0
        zero_to_one = lambda v: isinstance(v, float) and v >= 0 and v <= 1
        valid_number = lambda v: positive_int(v) or zero_to_one(v)

        for group, data in self.groups.items():
            min_in_groups[group] = 0
            for datum in data:
                count = datum.target_dist.total()
                min_in_groups[group] = min(count, min_in_groups[group])
                min_all_groups = min(count, min_all_groups)
        for group, data in self.groups.items():
            for datum in data:
                if datum.target_count == "min_in_groups": datum.target_count = min_in_groups[group]
                elif datum.target_count == "min_all_groups": datum.target_count = min_all_groups
                else: assert valid_number(datum.target_count), f"{datum.filename} had invalid target count {datum.target_count}"

    def downsample_to_target_count(self):
        for group, data in self.groups.items():
            for datum in data:
                datum.target_dist = datum.target_dist.to_count(datum.target_count)
    
    def create_samplers(self):
        for group, data in self.groups.items():
            for datum in data:
                datum.create_samplers()
    
    def write_samples(self):
        for group, data in self.groups.items():
            for datum in data:
                new_file = Path(datum.filename)
                parent = Path(datum.filename).parent
                new_filename = "downsample_" + Path(datum.filename).name
                print(new_filename)
                new_path = parent / new_filename
                print(new_path)
                datum.write_samples(str(new_path))

    def __str__(self):
        s = []
        for group, data in self.groups.items():
            s.append(f"Group {group}")
            for datum in data:
                s.append(str(datum))
        return "\n".join(s)



@click.command
@click.option("--groups", "-g", "--sample_groups", type = str, default = {})
@click.option("--downsample-to", "--to", type=DownsampleOption(), default = "min_all_groups")
@click.option("--outlier", multiple = True)
@click.option("--ignore-contig", "--ic", multiple = True)
@click.option("--ignore-cis", type = bool, default = False, show_default = True)
@click.option("--ignore-trans", type = bool, default = False, show_default = True)
@click.option("--cis-thresholds", "--thresholds", "--cis-strata", "--strata", type = str, default = default_cis_thresholds, show_default = True)
@click.option("--sep", default = "\t")
@click.option("--list-sep", default = ",")
def coverage(groups, downsample_to, outlier, ignore_contig, ignore_cis, ignore_trans, cis_thresholds, sep, list_sep):
    groups = load_groups(groups, sep, list_sep)
    space = TransCisThresholds("not pair.ur()", [10, 20, 50, 100])
    comparison = PairsComparison(space, groups, downsample_to)
    comparison.compute_counts() # Heavy processing - read through all lines in all files, also can probably be a separate function
    comparison.compute_central_distributions()
    comparison.downsample_to_target_probdist()
    comparison.compute_target_count()
    comparison.downsample_to_target_count()
    comparison.create_samplers()
    comparison.write_samples() # Heavy processing - write all sampled lines to all files

from collections import Counter

from numpy import searchsorted
from pathlib import Path

def pair_stats(pairs_file: PairsFile, output = None, columns = list[str], cis_strata = None):
    events = Counter()

    cols = ', '.join(columns)
    event_code = f"({cols})"
    event_code_compiled = compile(event_code, '<string>', 'eval')

    nat_partition = lambda cuts: sorted(list(set(cuts + [float('inf')]))) if cuts else None
    strata = nat_partition(cis_strata)

    if strata:
        get_stratum = lambda pair: strata[searchsorted(strata, pair.distance)] if pair.is_cis() else ""
    else:
        get_stratum = lambda pair: None

    for i, pair in enumerate(pairs_file):
        stratum = get_stratum(pair)
        event = eval(event_code_compiled)
        events[event] += 1
        if i > 10000:
            break
    
    stats_dict = defaultdict(list)
    for event, count in events.most_common():
        for col, row in zip(columns, event):
            stats_dict[col].append(str(row))
        stats_dict["count"].append(count)
    stats_df = pl.DataFrame(stats_dict)
    if output is None: print(stats_df)
    else: stats_df.write_csv(output, separator = "\t")
    

if __name__ == "__main__":
    pairs_file = PairsFile("data/240802_Arima_reseq_dedup.pairs.gz")
    stats(pairs_file, "stats.tsv", columns = ["pair.chr1", "pair.chr2", "stratum"], cis_strata = [10, 20, 50, 100])

"""
For stats, a space should define a way to map pair objects to a dict object
That dict gets 
"""

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

Architecture for setting up the run

Filename -> Pairs file + Histogram + Group
Set of groups ->
    -> Compute histograms
    -> Compute central histograms and assign as target histograms for files
    -> Adjust histogram to target profile
    -> Adjust histogram to specific size/fraction/sample size
    -> Downsample based on histogram

1. PairsSpace, {filename: group} -> {filename}
2. Extract 

"""

