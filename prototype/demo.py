
from hich.parse.pairs_file import PairsFile
from hich.stats.discrete_distribution import DiscreteDistribution
from hich.stats.pairs_rv import PairsRV
from hich.coverage.pairs_sampler import PairsSampler
from statistics import mean
import pypeln.process as pr
import pypeln.thread as pt
from pypeln.process import run

from pathlib import Path
import os
from dataclasses import dataclass, field
from typing import Callable
from collections import defaultdict
from itertools import chain
import polars as pl

@dataclass
class Collector:
    find: list
    extract: Callable
    collection: dict = field(default_factory = dict)

    def __call__(self, item):
        found = self.extract(item)
        if found in self.find:
            self.collection[found] = item
            if self.ready():
                yield self.items()

    def ready(self): return all([to_find in self.found() for to_find in self.find])
    
    def found(self): return list(self.collection.keys())

    def items(self): return list(self.collection.values())

@dataclass
class Grouper:
    groups: list[tuple]
    extract: Callable
    collection: dict = field(default_factory = dict)

    def __call__(self, item):
        ID = self.extract(item)
        for group in self.groups:
            if ID in group:
                self.collection.setdefault(group, dict())
                if ID not in self.collection[group]:
                    self.collection[group][ID] = item
                    if self.ready(group):
                        items = [item for item in self.collection[group].values()]
                        yield items
    
    def ready(self, group):
        found = self.collection[group].keys()
        return set(found) == set(group)


class View:
    def __call__(self, item):
        print(item)
        yield item

@dataclass
class Closure:
    closure: Callable

    def __call__(self, item):
        yield self.closure(item)

@dataclass
class ToCount:
    count: int

    def __call__(self, item):
        yield item.to_count(self.count)

@dataclass
class Flatten:
    def __call__(self, items):
        for item in items:
            yield item

@dataclass
class DownsamplePairsFile:
    input_path: str
    sampler: PairsSampler
    ID: int
    read_stats: str = ""
    write_stats: str = ""
    output_path: str = ""
    outlier: bool = False

    def __post_init__(self):
        if self.output_path == "":
            path = Path(self.input_path)
            parent = path.parent
            name = "downsampled_" + path.name
            self.output_path = str(parent / name) 

    def __str__(self):
        return f"ID: {self.ID}\n\tOutput: {self.output_path}\n\t{self.sampler.full}\n\t{self.sampler.target}\n"

@dataclass
class ComputePairsStats:
    max_records: int = float('inf')

    def __call__(self, item: DownsamplePairsFile):
        
        if Path(item.read_stats).expanduser().exists():
            read_stats = str(Path(item.read_stats).expanduser())
            df = pl.read_csv(read_stats, separator = "\t", has_header = True)
            item.sampler.full = item.sampler.rv.from_polars(df)
        else:
            pairs_file = PairsFile(item.input_path)
            item.sampler.count_events(pairs_file, self.max_records)
            if item.write_stats:
                df = item.sampler.rv.to_polars(item.sampler.full)
                df.write_csv(item.write_stats, separator = "\t", include_header = True)
        item.sampler.target = item.sampler.full.copy()
        yield item

class ToMeanMass:
    def __call__(self, items: list[DownsamplePairsFile]):
        targets = [item.sampler.target for item in items if not item.outlier]
        mean_mass = DiscreteDistribution.mean_mass(targets)
        for i in range(len(items)):
            items[i].sampler.target = items[i].sampler.target.downsample_to_probabilities(mean_mass)
        yield items

class ToMinGroupCount:
    def __call__(self, items: list[DownsamplePairsFile]):
        target_distributions = [item.sampler["target"] for item in items]
        min_count = min(target_distributions)
        for i, item in enumerate(items):
            items[i].sampler["target"] = items[i].sampler["target"].to_count(min_count)
        yield items

@dataclass
class WritePairsSample:
    max_records: int = float('inf')

    def __call__(self, item):
        pairs_file = PairsFile(item.input_path)
        item.sampler.write_sample(pairs_file, item.output_path, self.max_records)
        yield item

@dataclass
class ToCount:
    count: int

    def __call__(self, item: DownsamplePairsFile):
        item.sampler.target = item.sampler.target.to_count(self.count)
        yield item

def run_downsample_pipeline(files: list[DownsamplePairsFile],
               groups: list[tuple],
               mean_mass: bool = False,
               count: str = None,
               n_records: int = float('inf')):
    count = str(count)

    get_id = lambda item: item.ID
    compute_pairs_stats = ComputePairsStats(n_records)
    
    flatten = Flatten()
    
    to_group_min = ToMinGroupCount()
    view = View()

    write_sample = WritePairsSample(n_records)

    pipeline = (files | pr.flat_map(view) | pr.flat_map(compute_pairs_stats))

    if groups:
        to_sample_groups = Grouper(groups, lambda item: item.ID)
        pipeline = pipeline | pr.flat_map(to_sample_groups)
    
    if mean_mass:
        to_mean_mass = ToMeanMass()
        pipeline = pipeline | pr.flat_map(to_mean_mass)
    
    if count == "group_min":
        pipeline = pipeline | pr.flat_map(to_group_min)

    pipeline = pipeline | pr.flat_map(flatten)
    
    if count == "min":
        pipeline = pipeline | pr.flat_map(to_group_min)
    elif count.isnumeric():
        to_count = ToCount(int(count))
        pipeline = pipeline | pr.flat_map(to_count)
    elif count:
        try:
            pipeline = pipeline | pr.flat_map(float(to_count))
        except ValueError:
            pass
    
    pipeline = pipeline | pr.flat_map(write_sample) | pr.flat_map(view)

    run(pipeline)

wd = "~/Documents/hich/prototype/data/"
rv = PairsRV(["record.chr1", "record.chr2", "record.pair_type", "stratum"], [10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
file1 = DownsamplePairsFile(wd + "file1.pairs.gz", PairsSampler(rv = rv), ID = 1, read_stats = wd + "file1.stats", write_stats = wd + "file1.stats")
file2 = DownsamplePairsFile(wd + "file2.pairs.gz", PairsSampler(rv = rv), ID = 2, read_stats = wd + "file2.stats", write_stats = wd + "file2.stats")
file3 = DownsamplePairsFile(wd + "file3.pairs.gz", PairsSampler(rv = rv), ID = 3, read_stats = wd + "file3.stats", write_stats = wd + "file3.stats")
file4 = DownsamplePairsFile(wd + "file4.pairs.gz", PairsSampler(rv = rv), ID = 4, read_stats = wd + "file4.stats", write_stats = wd + "file4.stats")
files = [file1, file2, file3, file4]
groups = [(1, 2), (3, 4)]
run_downsample_pipeline(files, groups, True, 1000, 1000000)


