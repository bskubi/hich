import pypeln as pln
from pathlib import Path
from dataclasses import dataclass, field
from typing import Callable, Union
import os
from abc import ABC, abstractmethod
from hich.stats.discrete_distribution import DiscreteDistribution

class Nothing: ...

class PipelineStep(ABC):
    @abstractmethod
    def __call__(self): ...

def get_id(it): return it["id"]

@dataclass
class Collector(PipelineStep):
    collect_function: Callable = get_id
    to_collect: list = field(default_factory=list)
    collection: list = field(default_factory=dict)

    def __call__(self, item: dict) -> Union(Nothing, list[dict]):
        """Collect needed items and return collection when the last one is found

        Args:
            list (dict): A dict-like from which the ID can be extracted with
            the Collector's collect_function.

        Returns:
            Union: The list of collected items if the last one was just found,
            Nothing otherwise.
        """
        # Get ID and determine if the item is collectable.
        ID = self.get_id(item)
        item_is_collectable = ID in self.to_collect
        last_item_found = False
        if item_is_collectable:
            # If item is collectable, store it.
            self.collection[ID] = item

            # If the last item was just found, then set a flag to return the
            # collection of items.
            last_item_found = self.collected()

        # Return the collection of items if the last item was just found,
        # None otherwise.
        return self.collection if last_item_found else Nothing

    def collected(self) -> bool:
        """Return whether or not every ID to be collected has been collected

        Returns:
            bool: Whether or not all IDs to be collected have been collected
        """
        return all([ID in self.collection for ID in self.to_collect])

def get_distribution(items: Union(dict, list[dict])) -> list[dict]: ...


def downsample_to_central_distribution(items: Union(dict, list[dict]),
                                       center_function: Callable,
                                       ignore_key: str = "outlier"
                                       ) -> list[dict]:
    if not isinstance(items, list): items = [items]
    not_an_outlier = lambda item: not item[ignore_key]
    non_outliers = list(filter(not_an_outlier, items))
    central_distribution = center_function(non_outliers)
    downsample_item = lambda item: item.downsample_to_probabilities(central_distribution)
    downsampled = list(map(downsample_item, items))
    return downsampled

def downsample_to_count(items: Union(dict, list[dict]),
                        count: Union(str, int, float)):
    if not isinstance(items, list): items = [items]

    def item_count(item: dict) -> int:
        return sum(item["distribution"])

    min_count = min(map(item_count, items))
    target_count = min_count if count == "min" else count
    downsample_item = lambda item: item.to_count(target_count)
    downsampled = list(map(downsample_item, items))
    return downsampled