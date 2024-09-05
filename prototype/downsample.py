from collections import defaultdict
from dataclasses import dataclass, field
from hich.coverage.pairs_sampler import PairsSampler
from hich.parse.pairs_file import PairsFile
from hich.stats.discrete_distribution import DiscreteDistribution
from hich.stats.pairs_category import PairsCategorizer
from itertools import chain, filterfalse
from pathlib import Path
from pypeln.process import run
from statistics import mean
from typing import Callable, Generator, Tuple
import os
import polars as pl
import pypeln.process as pr
import pypeln.thread as pt

@dataclass
class Collector:
    """Collect required items based on identity and yield collection when ready."""
    collect: list
    identify: Callable
    collection: dict = field(default_factory = dict)

    def __call__(self, item: object) -> Generator:
        """Collect the required items based on their identity and yield the collection when ready

        Args:
            item (object): An item that may be collected based on its identity

        Yields:
            Generator: Yields a list of items when the last one is collected, or returns None otherwise.
        
        """

        # Use the identify function to get the identity of the item. If it should
        # be collected, then add it to the collection. If adding it completes
        # the collection, yield the collection. Otherwise return None.
        identity = self.identify(item)
        if identity in self.collect:
            self.collection[identity] = item
            if self.ready():
                yield self.items()

    def ready(self) -> bool:
        """Whether or not all items to be collected have been found"""
        return all([identity_to_collect in self.collected_identities() for identity_to_collect in self.collect])
    
    def collected_identities(self) -> list[object]:
        """Returns a list of the identities collected so far"""
        return list(self.collection.keys())

    def items(self):
        """Returns a list of the items collected so far"""
        return list(self.collection.values())

@dataclass
class Grouper:
    """Accumulate subgroups from a string of objects and emit subgroups when ready
    
    fields:
        groups (list[tuple]): A list of subgroup-defining tuples containing identifiers
                              for each subgroup item.
        identify (Callable): A function mapping an item to its identifier
        collection (dict):
            Key: a tuple of subgroup identifiers
            Value: A list of collected items in the order they were collected
    """
    groups: list[tuple]
    identify: Callable
    collection: dict = field(default_factory = dict)

    def __call__(self, item: object) -> Generator:
        """Accumulate subgroups and emit when ready.

        Subgroups are deemed ready to yield when at an instance of each
        subgroup item identifier has been collected. If the same ID is seen
        twice, all subgroups it belongs to will contain the most recent item
        associated with that ID.

        Args:
            item (object): An item whose identity can be extracted from the
                           assigned identify method. Will be stored with any
                           groups containing its ID.

        Yields:
            Generator: Yields a list of subgroup items when all items collected,
                        and returns None.
        """
        # Extract ID from item
        ID = self.identify(item)

        # Update the collection with the new ID and item
        groups_needing_ID = filter(lambda group: ID in group, self.groups)

        unready_groups_needing_ID = filterfalse(lambda group: self.ready(group),
                                                groups_needing_ID)

        map(self.add_item_to_group, unready_groups_needing_ID)

        # Yield any newly ready groups
        newly_ready_groups = filter(lambda group: self.ready(group),
                                    unready_groups_needing_ID)

        for group in newly_ready_groups:
            yield self.group_items(group)
    
    def add_item_to_group(self, group, ID, item):
        print(group, ID, item)
        self.collection.setdefault(group, dict())
        self.collection[group][ID] = item

    def group_items(self, group):
        return [item for item in self.collection[group].values()]
        

    def ready(self, subgroup_IDs: Tuple) -> bool:
        """Returns whether all members of a specific subgroup have been found.

        Args:
            group (_type_): _description_

        Returns:
            _type_: _description_
        """
        unique_collected_IDs = set(self.collection[subgroup_IDs].keys())
        unique_subgroup_IDs = set(subgroup_IDs)
        return unique_collected_IDs == unique_subgroup_IDs


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

grouper = Grouper([(1, 2), (3, 4)], lambda item: item["ID"])
next(grouper({"ID":1, "val": "one"}))

#

# wd = "~/Documents/hich/prototype/data/"
# rv = PairsRV(["record.chr1", "record.chr2", "record.pair_type", "stratum"], [10, 20, 50, 100, 200, 500, 1000, 2000, 5000])
# file1 = DownsamplePairsFile(wd + "file1.pairs.gz", PairsSampler(rv = rv), ID = 1, read_stats = wd + "file1.stats", write_stats = wd + "file1.stats")
# file2 = DownsamplePairsFile(wd + "file2.pairs.gz", PairsSampler(rv = rv), ID = 2, read_stats = wd + "file2.stats", write_stats = wd + "file2.stats")
# file3 = DownsamplePairsFile(wd + "file3.pairs.gz", PairsSampler(rv = rv), ID = 3, read_stats = wd + "file3.stats", write_stats = wd + "file3.stats")
# file4 = DownsamplePairsFile(wd + "file4.pairs.gz", PairsSampler(rv = rv), ID = 4, read_stats = wd + "file4.stats", write_stats = wd + "file4.stats")
# files = [file1, file2, file3, file4]
# groups = [(1, 2), (3, 4)]
# run_downsample_pipeline(files, groups, True, 1000, 1000000)


