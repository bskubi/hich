from .df_required_cols import DFRequiredCols
from .frag_index import FragIndex
import polars as pl

class BedpeIDEnds(DFRequiredCols):
    required_cols = ["pairID", "end", "chrom", "pos"]
    
    def __init__(self, df, prev_names = {}):
        super().__init__(df)
        self.prev_names = prev_names

    def reformat(self, fmt):
        from .bedpe_id_pairs import BedpeIDPairs
        BedpeIDEnds.conversions = [BedpeIDPairs]
        self.assert_conversion_allowed(fmt)

        if fmt == BedpeIDPairs:
            ends = {}
            for end in ["1", "2"]:
                frag_cols = {
                    "rfrag": f"rfrag{end}",
                    "rfrag_start": f"rfrag_start{end}",
                    "rfrag_end": f"rfrag_end{end}"
                }
                end_cols = self.prev_names[end]
                new_cols = {**frag_cols, **end_cols}
                
                df = self.df.filter(pl.col("end") == end) \
                            .drop("end") \
                            .rename(new_cols)

                ends[end] = df


            df = ends["1"].join(ends["2"],
                                on="pairID",
                                suffix = "__dropme__",
                                coalesce=True)
            drop_cols = [col for col in df.columns if col.endswith("__dropme__")]
            df = df.drop(drop_cols) \
                   .sort(by=["chrom1", "chrom2", "pos1", "pos2"])

            return BedpeIDPairs(df)

    def fragtag(self, frag_index):
        chrom_partition = self.df.partition_by(['chrom'], as_dict = True)
        for chrom, chrom_df in chrom_partition.items():
            chrom_partition[chrom] = self.get_frag_columns(frag_index,
                                                           chrom,
                                                           chrom_df)
        self.df = pl.concat(chrom_partition.values())
        return self

    def get_frag_columns(self,
                         frag_index,
                         chrom,
                         chrom_df,
                         frag_colnames = ["rfrag", "rfrag_start", "rfrag_end"]):
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

        frag_cols = [rfrag, rfrag_start, rfrag_end]

        for i, frag_col in enumerate(frag_cols):
            frag_cols[i] = frag_col.cast(pl.Int64).alias(frag_colnames[i])

        return chrom_df.with_columns(frag_cols)