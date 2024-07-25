import click
import numpy as np
import polars as pl
import time
import warnings
import sys
from .pairs_parser import PairsParser
from .frag_index import FragIndex
from .bedpe_pairs import BedpePairs
from .samheader_fragtag import SamheaderFragtag

def tag_batch(pairs_batch_df, frag_index):
    pairs_batch_df = tag_with_pair_id(pairs_batch_df)
    ends = ends_format(pairs_batch_df)
    chroms_dict = partition_by_chrom(ends)

    for chrom, chrom_df in chroms_dict.items():
        frag_columns = format_frag_columns(frag_index,
                                            chrom,
                                            chrom_df)
        chroms_dict[chrom] = chroms_dict[chrom].with_columns(frag_columns)
    
    return pairs_format(pairs_batch_df, chroms_dict)

def tag_restriction_fragments(frags_filename: str,
                              input_pairs_filename: str,
                              output_pairs_filename: str,
                              batch_size: int = 1000000):
    
    frag_index = FragIndex(frags_filename)

    pairs_parser = PairsParser(input_pairs_filename)

    for df in pairs_parser.batch_iter(batch_size):  
        df = BedpePairs(df).fragtag(frag_index)

        pairs_parser.write_append(output_pairs_filename,
                                  df,
                                  header_end = SamheaderFragtag())


    pairs_parser.close()

@click.command()
@click.option("--batch_size", default = 1000000)
@click.argument("fragfile")
@click.argument("out_pairs")
@click.argument("in_pairs")
def fragtag(batch_size, fragfile, out_pairs, in_pairs):
    tag_restriction_fragments(fragfile,
                            in_pairs,
                            out_pairs,
                            batch_size)

if __name__ == "__main__":
    fragtag()

# python fragtag /home/benjamin/Documents/hich2/work/8f/32dfb1ed645250689bf66d3c7664bb/hg38_Arima.bed result.pairs.gz /home/benjamin/Documents/hich2/work/8f/32dfb1ed645250689bf66d3c7664bb/Mock_2_1.pairs.gz