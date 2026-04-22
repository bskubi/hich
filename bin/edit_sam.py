#!/usr/bin/env python

import pysam
import sys
from hich.core.namespace import namespace_from_path
import click
from typing import Protocol

class SAM_Editor(Protocol):
    def edit_seg(self, seg: pysam.AlignedSegment) -> pysam.AlignedSegment: ...

@click.command
@click.option("--edit-seg-func", default="edit_seg")
@click.argument("config")
def sam_editor(edit_seg_func, config):
    ns = namespace_from_path(config, "edit_sam", SAM_Editor)
    decompression_threads = getattr(ns, "decompression_threads", 3)
    compression_threads = getattr(ns, "compression_threads", 8)
    sam_input = getattr(ns, "sam_input", "-")
    sam_output = getattr(ns, "sam_output", "-")
    write_mode = getattr(ns, "write_mode", "w")
    
    edit_seg = getattr(ns, edit_seg_func, lambda x: x)
    sam_input = pysam.AlignmentFile(sam_input, threads=decompression_threads)
    sam_output = pysam.AlignmentFile(sam_output, template=sam_input, mode=write_mode, threads=compression_threads)
    

    for seg in sam_input:
        sam_output.write(edit_seg(seg))

if __name__ == "__main__":
    sam_editor()