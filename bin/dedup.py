#!/usr/bin/env python

import subprocess
import click
"""
1. Edit SAM segments. Required because samtools markdup only accepts
one barcode tag, but some datasets have unique cell identifier over multiple
tags or located in the read ID.
2. Sort by name (required for fixmates)
3. Fix mates (required for markdup)
4. Sort by coordinate (required for markdup)
5. Mark duplicates and output duplicate statistics.
6. Sort by name (required for pairtools parse, no impact on methylation calling)
"""
@click.command
@click.option("--barcode-tag", type=str)
@click.option("--n-procs", type=int, default=1)
@click.option("--config-edit-sam", type=str, default="")
@click.option("--stats-output", type=str, default="stats.json")
@click.option("--sam-output", type=str, default="-")
def dedup(barcode_tag, n_procs, config_edit_sam, stats_output, sam_output):
    perf = f"-O BAM -u -@ {n_procs}"
    barcode_tag_arg = f"--barcode-tag {barcode_tag}" if barcode_tag else ""
    stats_output_arg = f"-f {stats_output}" if stats_output else ""
    sam_output_arg = f"> {sam_output}" if sam_output else ""

    cmd = " | ".join([
        f"edit_sam.py --edit-seg-func edit_seg1 {config_edit_sam}",
        f"samtools sort {perf} -n",
        f"samtools fixmate {perf} -m - - ",
        f"samtools sort {perf}",
        f"samtools markdup {perf} {stats_output_arg} {barcode_tag_arg} --json --duplicate-count -S --include-fails - -",
        f"samtools sort {perf} -n",
        f"edit_sam.py --edit-seg-func edit_seg2 {config_edit_sam} {sam_output_arg}"
    ])
    print(cmd)

    subprocess.run(cmd, shell=True)

if __name__ == "__main__":
    dedup()