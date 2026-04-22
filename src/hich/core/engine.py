"""
Batch processing engine

Read in chunks from data
Parallel process, return results
Compute stats on results
Store result in TileDB
Store stats
"""
from dataclasses import dataclass, field
from typing import Iterator
import polars as pl
from typing import TypeVar, Generic, Protocol, Optional, List, Callable, Dict
from .statistics import Statistic, StatisticType, StatisticAccumulator
from pathlib import Path
from .namespace import namespace_from_path
from typing import Any, Self
from abc import ABC, abstractmethod
import tiledb as tdb
from . import tiledb_utils as tdb_utils
import signal
from loky import ProcessPoolExecutor
import multiprocessing
import os
import time
from collections import defaultdict

class Batch:
    ...

@dataclass
class Result:
    observations: pl.DataFrame


class EngineNamespace(Protocol):
    def observations_schema(self) -> tdb.ArraySchema:
        """
        Get the schema for the observations
        """
        ...

    def statistics(self) -> Optional[Dict[str, StatisticType]]:
        """
        Get the schema for the statistics to be computed over the observations
        """
        ...

    def validate_observations(observations: pl.DataFrame) -> pl.DataFrame:
        """
        Validate output observations prior to writing to disk or computing statistics 
        """
        ...

    def write_statistics(self, name: str, df: pl.DataFrame):
        ...

def default_write_statistics(name: str, df: pl.DataFrame):
    df.write_csv(f"{name}.tsv", separator="\t")
    df.write_parquet(f"{name}.parquet")


T = TypeVar("T", bound=EngineNamespace)

@dataclass
class Engine(ABC, Generic[T]):
    namespace: T | None = None
    statistics_dir: Path = field(default_factory=lambda: Path("./"))
    statistic_types: Dict[str, StatisticType] | None = None
    statistic_accumulators: Dict[str, StatisticAccumulator] | None = None
    
    validate_observations: Callable[[pl.DataFrame], pl.DataFrame] = None
    write_statistics: Callable[[str, pl.DataFrame], None] | None = None
    
    observations_schema: tdb.ArraySchema | None = None
    max_workers: Optional[int] = None
    max_queue_depth: Optional[int] = None

    timing: Dict[str, float] = field(default_factory=lambda: defaultdict(float))

    def load_config_py(self, config_py_path: str | Path) -> Self:
        ns = namespace_from_path(config_py_path, "engine_ns", EngineNamespace)
        self.namespace = ns
        self.statistic_types = ns.statistics() or {}
        self.statistic_accumulators = {name: stat.accumulator(stat) for name, stat in self.statistic_types.items()}
        self.observations_schema = ns.observations_schema()
        self.validate_observations = ns.validate_observations if hasattr(ns, "validate_observations") else lambda x: x
        if hasattr(ns, "write_statistics"):
            self.write_statistics = ns.write_statistics
        else:
            self.write_statistics = default_write_statistics
        
        return self
    
    def run(self, input_stream: Any, observations_path: Path, append_observations: bool=False):
        """
        Run pipeline to process batches, 
        """
        obs_array = self._open_observations_array(observations_path, append_observations)
        
        for result in self._process_batches(input_stream):   
            start_time = time.time()
            self._handle_observations_result(result, obs_array)
            self.timing["handle_observations_result"] += time.time() - start_time
            start_time = time.time()
            self._handle_statistics_result(result)
            self.timing["handle_statistics_result"] += time.time() - start_time
        start_time = time.time()
        self.finish_writes(obs_array)
        self.timing["finish_writes"] = time.time() - start_time

    @property
    def ns(self) -> T:
        return self.namespace

    @staticmethod
    @abstractmethod
    def worker(self, *args, **kwargs) -> Result:
        ...

    @abstractmethod
    def batches(self, input_stream: Any) -> Iterator[dict]:
        """
        Batch input_stream and yield dicts containing arguments to self.worker
        """
        ...

    def _process_batches(self, input_stream: Any) -> Iterator[Result]:
        batches = self.batches(input_stream)
        max_queue_depth = self.max_queue_depth or os.cpu_count() * 2

        with ProcessPoolExecutor(
            max_workers=self.max_workers, 
            initializer=init_worker
        ) as ppe:

            # Collect a set of futures 
            futures = []
            
            def submit_next_batch() -> bool:
                
                try:
                    start_time = time.time()
                    batch = next(batches)
                    self.timing["get_next_batch"] += time.time() - start_time
                except StopIteration:
                    return False

                # Submit another process to be run when a CPU is available.
                future = ppe.submit(self.worker, **batch)
                futures.append(future)
                return True

            # Pre-fill execution queue
            for _ in range(max_queue_depth):
                if not submit_next_batch():
                    break
                    
            # Drain and replenish queue sequentially
            
            while futures:
                start_time = time.time()
                result: Result = futures.pop(0).result()
                result.observations = self.validate_observations(result.observations)
                yield result
                self.timing["futures.pop(0).result()"] += time.time() - start_time
                submit_next_batch()

    def _open_observations_array(self, path: Path, append: bool) -> tdb.SparseArray:
        tdb_utils.create(uri=str(path), schema=self.observations_schema, append=append)
        obs_array = tdb.SparseArray(str(path), "w")
        return obs_array

    def _handle_statistics_result(self, result: Result):
        for stat_accum in self.statistic_accumulators.values():
            stat_accum.update(result.observations)

    def _handle_observations_result(self, result: Result, array: tdb.SparseArray):
        tdb_utils.write_pl_to_tiledb(result.observations, array)

    def finish_writes(self, array: tdb.Array):
        array.close()
        stats_dir = Path(self.statistics_dir)
        stats_dir.mkdir(parents=True, exist_ok=True)
        for name, stat in self.statistic_accumulators.items():
            self.write_statistics(str(stats_dir / name), stat.statistic.df)



def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)