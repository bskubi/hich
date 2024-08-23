from hich.stats.discrete_distribution import DiscreteDistribution
from hich.coverage.pairs_sampler import PairsSampler
from statistics import mean
import pypeln as pln
from pypeln.process import run
from hich.parse.pairs_file import PairsFile
from pathlib import Path
import os
from dataclasses import dataclass, field
from typing import Callable
from collections import defaultdict
from itertools import chain


def mp(function, *args, **kwargs):
    return pln.process.map(function, *args, workers=os.cpu_count(), **kwargs)

@dataclass
class DownsampleStep:
    ids: list
    size: int = None
    outliers: list = field(default_factory=list)
    center: Callable = None
    distributions: defaultdict = field(default_factory=lambda: defaultdict(lambda: None))
    target: DiscreteDistribution = None

    def __call__(self, data: DiscreteDistribution) -> dict[object, DiscreteDistribution]:
        if isinstance(data, DownsampleStep):
            for distribution in data.distributions:
                self(distribution)

        ID = data["id"]
        if ID in self.ids:
            self.distributions[ID] = DiscreteDistribution(data["distribution"])
            if self.collected():
                self.downsample()
                return self.result()

    def downsample(self):
        non_outliers = [dist for ID, dist in self.distributions.items() if ID not in self.outliers]
        if self.center:
            self.target = self.center(non_outliers)
            for ID, distribution in self.distributions.items():
                self.distributions[ID] = distribution.downsample_to_probabilities(self.target)

        self.size = min(self.distributions.values()).total() if self.size == "min" else self.size
        if self.size:
            for ID, distribution in self.distributions.items():
                self.distributions[ID] = distribution.to_count(self.size)
        
    def collected(self):
        return all([self.distributions[ID] is not None for ID in self.ids])
    
    def result(self):
        return [{"id": ID, "distribution": distribution} for ID, distribution in self.distributions.items()]

def collect(pipeline):
    return [result for result in list(pipeline) if result is not None]

def compute_downsample_targets(filename_subgroups: list[list[str]],
               center: bool,
               size: str,
               conjuncts = list[str],
               cis_strata = None):

    def float_or_none(f):
        f = str(f)
        try: return int(f) if f.isnumeric() else float(f)
        except ValueError: return None
    
    sampler = PairsSampler(conjuncts, cis_strata)
    get_stats = StatsCaller(sampler, 10000)
    filenames = list(chain.from_iterable(filename_subgroups))
    stats = list(filenames | mp(get_stats))
    center = DiscreteDistribution.mean_mass if center else None
    size1 = "min" if size == "min_subgroup" else None
    size2 = "min" if size == "min_all" and size != "min_subgroup" else float_or_none(size)
    
    step1 = [DownsampleStep(subgroup, center = center, size = size1) for subgroup in filename_subgroups]
    step2 = DownsampleStep(filenames, size = size2)
    
    results = []

    for subgroup_target in step1:
        result = collect(stats | mp(subgroup_target))[0]
        results.extend(result)

    results = collect(results | mp(step2))[0]

    return results

# There should be a tool that converts records into events that you load with cis_strata and conjuncts.
# This way you could create the object once, and use it both for computing stats and for sampling
# It could be used to compute stats directory (counting up) and also for sampling (counting down)

filenames = [f"data/file{i}.pairs.gz" for i in range(1, 5)]
sample_groups = [[filenames[0], filenames[1]], [filenames[2], filenames[3]]]
print(compute_downsample_targets(sample_groups, False, .5, ["pair.chr1", "pair.chr2", "pair.pair_type", "stratum"], [100, 200, 500, 1000, 2000, 5000]))

"""
Architecture:

DiscreteDistribution
SelectionSampler

PairsRV: pair -> event
PairsSampler(SelectionSampler): compute statistics over and sample pairs files

DownsampleStep
DownsamplePipeline: PairsFiles + parameters -> PairsSampler for each file
"""