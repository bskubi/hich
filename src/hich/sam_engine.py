import os
import sys
import subprocess
from pathlib import Path
from dataclasses import dataclass
from typing import Iterator, Tuple, IO, Any
from hich.core import Engine
from abc import ABC, abstractmethod

def fast_stream(path, size, cpus):
    # Isolated header capture
    h = subprocess.check_output(["samtools", "view", "-H", "-@", str(cpus), str(path)])
    # Alignment-only stream
    p = subprocess.Popen(["samtools", "view", "-@", str(cpus), str(path)], stdout=subprocess.PIPE)
    
    r = b""
    while c := p.stdout.read(size):
        # Optimization: Assume records < size; avoid conditional accumulation
        i = c.rfind(b'\n') + 1
        yield h, r, c[:i]
        r = c[i:]
        
    if r: yield h, r, b""
    p.kill()

@dataclass
class SAM_Engine(Engine):
    batch_size_bytes: int = 512 * 1024 * 1024 
    samtools_view_cpus: int = 16

    def batches(self, input_stream: Path | str) -> Iterator[dict]:
        if Path(input_stream).exists():
            for h, b, t in fast_stream(input_stream, self.batch_size_bytes, self.samtools_view_cpus):
                yield self.SAM_batches({"sam_header": h, "sam_body": b, "sam_tail": t })
        else:
            raise NotImplementedError("input_stream must be a path")

    def SAM_batches(self, sam_stream: dict) -> dict:
        return sam_stream

    @classmethod
    @abstractmethod
    def SAM_worker(cls, sam_stream: bytes, *args, **kwargs) -> Any:
        ...
    
    @classmethod
    def worker(cls, sam_header: bytes, sam_body: bytes, sam_tail: bytes, *args, **kwargs) -> Any:
        return cls.SAM_worker(sam_stream = sam_header + sam_body + sam_tail, *args, **kwargs)