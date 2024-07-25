from .df_required_cols import DFRequiredCols
import polars as pl

class BedpeIDPairs(DFRequiredCols):
    required_cols = ["pairID"]
    
    def __init__(self, df):
        super().__init__(df)
    
    def reformat(self, fmt):
        from .bedpe_id_ends import BedpeIDEnds
        from .bedpe_pairs import BedpePairs
        BedpeIDPairs.conversions = [BedpeIDEnds, BedpePairs]
        self.assert_conversion_allowed(fmt)

        if fmt == BedpeIDEnds:
            ends = []
            end_names = {}
            for end, other in [("1", "2"), ("2", "1")]:
                drop = [col
                        for col in self.df.columns
                        if col.endswith(other)]
                
                rename = {col: col.strip(end)
                          for col in self.df.columns
                          if col.endswith(end)}
                
                end_names[end] = {val: key for key, val in rename.items()}
                df = self.df.drop(drop) \
                            .rename(rename) \
                            .with_columns(pl.lit(end).alias("end"))
                ends.append(df)

            df = pl.concat(ends)
            
            return BedpeIDEnds(df, end_names)


        elif fmt == BedpePairs:
            from .bedpe_pairs import BedpePairs
            df = self.df.drop("pairID")
            return BedpePairs(df)