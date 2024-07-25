import polars as pl
import numpy as np

def tag_with_pair_id(pairs, colname = "pairID"):
    # We will be sorting each position separately, so to keep track of which pair
    # each end originally belonged to we assign it a pair_id
    pairs_in_batch = len(pairs)
    pair_ids = range(pairs_in_batch)
    pairs = pairs.with_columns(pl.Series(pair_ids).alias("pairID"))
    return pairs

def frag_default_cols():
    # Pairtools uses rfrag1/2, rfrag_start1/2, and rfrag_end1/2 as column names
    # And -1, 0, 0 respectively as fillers for unmapped reads, so we create default
    # columns to be renamed with the 1/2 suffix later.
    frag_cols = ['rfrag', 'rfrag_start', 'rfrag_end']
    unmapped_default_entries = [-1, 0, 0]
    frag_default_cols = [pl.lit(unmapped_default_entry) \
                           .alias(colname) \
                           .cast(pl.Int64)
                         for colname, unmapped_default_entry
                         in zip(frag_cols, unmapped_default_entries)]

    return frag_default_cols

def ends_format(pairs):
    # For each end, sort by chrom and pos, then break into a per-chrom dict
    # Final columns are: pair_id, chrom, pos, rfrag, rfrag_start, rfrag_end
    ends = []
    
    for end in [1,2]:
        # Select pair_id, chrom, and pos columns
        select_col_names = ['pairID', f'chrom{end}', f'pos{end}']
        
        # Convert chromosome and position with end identifier numbers
        # to generic 'chrom' and 'pos'
        new_col_names = {f'chrom{end}': 'chrom', f'pos{end}':'pos'}
        
        # Marker for the pair end (1 or 2)
        end_col = pl.repeat(end, len(pairs)).alias("end")

        # Select and rename columns, then add default frag columns and end identifier
        end_df = pairs.select(select_col_names) \
                    .rename(new_col_names) \
                    .with_columns(frag_default_cols()) \
                    .with_columns(end_col)
        ends.append(end_df)

    # Combine into single dataframe and sort by chrom and position
    ends = pl.concat(ends)
    return ends

def partition_by_chrom(ends):
    # Partition by chromosome into a dictionary {(chrom,): DataFrame}
    # With pair_id, pos, rfrag, rfrag_start, rfrag_end, and end columns
    return ends.partition_by(['chrom'],
                            as_dict = True,
                            include_key = False)


def format_frag_columns(frag_index, chrom, chrom_df):
    # Find restriction fragment indices for each read
    positions = chrom_df['pos'].to_list()

    if chrom not in frag_index:
        count = len(chrom_df)
        rfrag = pl.repeat(-1, count)
        rfrag_start = pl.repeat(0, count)
        rfrag_end = pl.repeat(0, count)
    else:
        # Get the indices where the positions sort to in the restriction fragments
        # for the chromosome -- this is the main event that identifies the restriction fragments
        rfrag_indices = frag_index.search(chrom, positions)

        # Create column for the restriction fragment index and for the start and end positions
        rfrag = pl.Series(rfrag_indices)

        # Create column for start position of the restriction fragment
        # Adjust exact positioning to match pairtools output
        # (Set 0 to -1, then add 1 to all)
        rfrag_start = frag_index.starts(chrom) \
                        .gather(rfrag_indices) \
                        .replace(0, -1) \
                        + 1

        # Create column for end position of the restriction fragment
        # Adjust exact position to match pairtools output (add 1)
        rfrag_end = frag_index.ends(chrom) \
                    .gather(rfrag_indices) \
                    + 1
    rfrag = rfrag.cast(pl.Int64).alias('rfrag')
    rfrag_start = rfrag_start.cast(pl.Int64).alias('rfrag_start')
    rfrag_end = rfrag_end.cast(pl.Int64).alias('rfrag_end')

    # Replace the default rfrag columns with the new info
    return [rfrag, rfrag_start, rfrag_end]

def pairs_format(df, chroms_dict):
    # Rejoin fragments with main dataframe
    ends = {1:[], 2:[]}
    new_colnames = {}
    for chrom, chrom_df in chroms_dict.items():
        ends_dict = chrom_df.partition_by(['end'],
                                as_dict = True,
                                include_key = False)

        for end_key, end_dict in ends_dict.items():
            # Extract integer as the partition_by method of Polars stores the end
            # as a single-item tuple
            end = end_key[0]
            ends[end].append(end_dict)
        
            new_colnames[end] = {'rfrag':       f'rfrag{end}', 
                                 'rfrag_start': f'rfrag_start{end}',
                                 'rfrag_end':   f'rfrag_end{end}'}

    
    ends = [pl.concat(ends[end]) \
                .rename(new_colnames[end]) \
                .drop('pos')
            for end
            in [1,2]]

    for end in ends:
        df = df.join(end, on = 'pairID', how = 'left')

    df = df.drop('pairID')
    return df

