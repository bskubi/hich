import click
import numpy as np
import polars as pl
import time
import polars as pl
import smart_open

class PairsParser:
    def __init__(self, filename):
        self.filename = filename
        self.write_file = None
        
    def columns_row(self):
        """Returns the 0-indexed offset of the final header row with the column names

        row_index: 0-indexed offset of the final header row (header rows start with '#')
        columns_dict: {colname:index} dict of which columns correspond to which column names

        If a row starting with '#columns:' is present it is assumed to be the final header row.
        If no '#columns:' row is present, returns the last header row with the minimal seven
        column names given by the .pairs specification.
        
        Returns: (row_index, columns_dict)
        """

        default_columns = ['readID', 'chrom1', 'pos1', 'chrom2', 'pos2', 'strand1', 'strand2']
        row = 0

        with smart_open.open(self.filename, "rt", encoding='utf-8') as file:
            for line in file:
                if line.startswith('#columns:'):
                    column_names = line.split()[1:]
                    return (row, column_names)
                elif line.startswith('#'):
                    row += 1
                else:
                    return (row-1, default_columns)
        return (row - 1, default_columns)

    def header(self):
        """Returns list of all lines starting with '#' until but not including the first line not starting with '#'
        """
        header = []
        with smart_open.open(self.filename, "rt", encoding="utf-8") as file:
            for line in file:
                if line.startswith('#'):
                    header.append(line)
                else:
                    break
        return header

    def batch_iter(self, n_rows):
        def read_chunk(skip_rows, columns):
            return pl.read_csv(self.filename,
                             skip_rows = skip_rows,
                             n_rows = n_rows,
                             raise_if_empty = False,
                             separator = '\t',
                             has_header = False,
                             new_columns = columns,
                             dtypes = {"readID":pl.String,
                                       "chrom1":pl.String,
                                       "chrom2":pl.String})
            
        columns_row, columns = self.columns_row()

        skip_rows = columns_row + 1
        with open(self.filename) as file:
            df = read_chunk(skip_rows, columns)
            skip_rows += len(df)
            yield df
            while len(df) == n_rows:
                df = read_chunk(skip_rows, columns)
                skip_rows += len(df)
                if not df.is_empty():
                    yield df

    def write_append(self, df, filename):
        import smart_open

        if self.write_file is None:
            self.write_file = smart_open.open(filename, "w")
            header_no_columns = "".join(self.header()[:-1])
            columns = " ".join(["#columns:"] + df.columns) + "\n"
            header = header_no_columns + columns
            self.write_file.write(header)

        df.write_csv(self.write_file,
                     include_header = False,
                     separator = '\t')

    def close(self):
        if self.write_file is not None:
            self.write_file.close()
            self.write_file = None

def tag_restriction_fragments(frags_filename: str,
                              input_pairs_filename: str,
                              output_pairs_filename: str,
                              batch_size: int = 1000000):

    ###################################################
    # README

    # The key idea here is to identify restriction fragments
    # by sorting the positions of each pair end separately,
    # as well as the position of each RE fragment end,
    # then using np.searchsorted() to find where in the list
    # of sorted RE fragment ends each position would fall.
    
    # Read zero-indexed .bed format restriction fragments file
    # Sort by chrom and start position
    # Then split into per-chromosome {(chrom,): DataFrame (start end)} dicts
    # Where the DataFrame includes only the start and end positions

    frags = pl.read_csv(frags_filename,
                        separator = '\t',
                        has_header = False,
                        new_columns = ['chrom', 'start', 'end']) \
               .partition_by(['chrom'],
                             as_dict = True,
                             include_key = False)


    # Cache the end positions in a list for repeated use when searching for restriction fragments
    frag_ends_dict = {chrom:
                      sorted(frags[chrom]['end'].to_list())
                      for chrom
                      in frags}

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

    # Parser to manage extracting batches of pairs from .pairs format
    pp = PairsParser(input_pairs_filename)

    # Iteratively extract batches of reads
    for df in pp.batch_iter(batch_size):  
        
        # We will be sorting each position separately, so to keep track of which pair
        # each end originally belonged to we assign it a pair_id
        pairs_in_batch = len(df)
        pair_ids = range(pairs_in_batch)
        df = df.with_columns(pl.Series(pair_ids).alias("pair_id"))

        # For each end, sort by chrom and pos, then break into a per-chrom dict
        # Final columns are: pair_id, chrom, pos, rfrag, rfrag_start, rfrag_end
        ends = []
        for end in [1,2]:
            # Select pair_id, chrom, and pos columns
            select_col_names = ['pair_id', f'chrom{end}', f'pos{end}']
            
            # Remove pair end identifier from column
            new_col_names = {f'chrom{end}': 'chrom', f'pos{end}':'pos'}
            
            # Marker for the pair end (1 or 2)
            end_col = pl.repeat(end, pairs_in_batch).alias("end")

            # Select and rename columns, then add default frag columns and end identifier
            end_df = df.select(select_col_names) \
                       .rename(new_col_names) \
                       .with_columns(frag_default_cols) \
                       .with_columns(end_col)
            ends.append(end_df)

        # Combine into single dataframe and sort by chrom and position
        ends = pl.concat(ends)

        # Partition by chromosome into a dictionary {(chrom,): DataFrame}
        # With pair_id, pos, rfrag, rfrag_start, rfrag_end, and end columns
        chroms_dict = ends.partition_by(['chrom'],
                                        as_dict = True,
                                        include_key = False)
    
        for chrom, chrom_df in chroms_dict.items():
            
            # Make sure we have restriction fragments for the current chromosome
            if chrom in frag_ends_dict:
                
                
                # Find restriction fragment indices for each read
                positions = chrom_df['pos'].to_list()

                # Get the indices where the positions sort to in the restriction fragments
                # for the chromosome -- this is the main event that identifies the restriction fragments
                rfrag_indices = np.searchsorted(frag_ends_dict[chrom], positions)
                
                if chrom == ('chrM',):
                    print(len(frags[chrom]['start']))
                    max_index = rfrag_indices.tolist().index(max(rfrag_indices))
                    print(positions[max_index])
                    print(frag_ends_dict[chrom])


                # Create column for the restriction fragment index and for the start and end positions
                rfrag = pl.Series(rfrag_indices) \
                          .cast(pl.Int64) \
                          .alias('rfrag')

                # Create column for start position of the restriction fragment
                # Adjust exact positioning to match pairtools output
                # (Set 0 to -1, then add 1 to all)
                rfrag_start = frags[chrom]['start'] \
                              .gather(rfrag_indices) \
                              .cast(pl.Int64) \
                              .alias('rfrag_start') \
                              .replace(0, -1) \
                              + 1

                # Create column for end position of the restriction fragment
                # Adjust exact position to match pairtools output (add 1)
                rfrag_end = frags[chrom]['end'] \
                            .gather(rfrag_indices) \
                            .cast(pl.Int64) \
                            .alias('rfrag_end') \
                            + 1

                # Replace the default rfrag columns with the new info
                frag_columns = [rfrag, rfrag_start, rfrag_end]
                chroms_dict[chrom] = chroms_dict[chrom].with_columns(frag_columns)
        
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
            df = df.join(end, on = 'pair_id', how = 'left')
            
        df = df.drop('pair_id', 'end')
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