import click
import numpy as np
import polars as pl
import time
import warnings
from .pairsParser import PairsParser
from .pairsTagger import *

def tag_restriction_fragments(frags_filename: str,
                              input_pairs_filename: str,
                              output_pairs_filename: str,
                              batch_size: int = 1000000):
    
    frags, frag_ends_dict = readFrags(frags_filename)

    pp = PairsParser(input_pairs_filename)

    for df in pp.batch_iter(batch_size):  
        df = tagWithPairID(df)
        ends = columnsPairIDChromPosEndFrags(df)
        chroms_dict = partitionEndsByChrom(ends)

        for chrom, chrom_df in chroms_dict.items():
            if chrom in frag_ends_dict:
                chroms_dict = labelFragments(frags, frag_ends_dict, chrom, chroms_dict, chrom_df)
        
        df = toPairsFormatDF(df, chroms_dict, chrom_df)

        pp.write_append(df, output_pairs_filename)
    pp.close()

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