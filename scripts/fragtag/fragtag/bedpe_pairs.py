import polars as pl

class BedpePairs():
    def __init__(self, df):
        self.df = df
    
    def fragtag(self, frag_index):
        pair_id = "__pairID__"
        pair_id_df = self.df.with_row_index(pair_id)
        

        for end in ["1", "2"]:
            chrom_col = f"chrom{end}"
            pos_col = f"pos{end}"
            frag_colnames = [f"rfrag{end}", f"rfrag_start{end}", f"rfrag_end{end}"]

            end_chroms = pair_id_df.select(pair_id, chrom_col, pos_col) \
                                   .partition_by([chrom_col], as_dict = True)

            for chrom, chrom_df in end_chroms.items():
                chrom_df = chrom_df.sort(by = [pos_col])
                positions = chrom_df[pos_col].to_list()
                frag_cols = BedpePairs.frag_columns(frag_index,
                                                    chrom,
                                                    positions,
                                                    frag_colnames)
                chrom_df = chrom_df.with_columns(frag_cols)
                chrom_df = chrom_df.drop([chrom_col, pos_col])
                frag_cols = chrom_df.sort(by = [pair_id]) \
                                    .select(frag_colnames)                              
                end_chroms[chrom] = end_chroms[chrom].with_columns(frag_cols)
            
            end_frags = pl.concat(end_chroms.values()) \
                          .drop([chrom_col, pos_col])

            pair_id_df = pair_id_df.join(end_frags,
                                         on = [pair_id],
                                         how = "inner",
                                         coalesce = True)
        pair_id_df = pair_id_df.sort(by = [pair_id])
        return pair_id_df.drop(pair_id)

    @classmethod
    def frag_columns(self,
                         frag_index,
                         chrom,
                         positions,
                         frag_colnames):
        if chrom not in frag_index:
            count = len(positions)
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
        frag_cols = [rfrag, rfrag_start, rfrag_end]
        for i, frag_col in enumerate(frag_cols):
            frag_cols[i] = frag_col.cast(pl.Int64).alias(frag_colnames[i])

        return frag_cols