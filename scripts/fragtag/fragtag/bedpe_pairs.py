from .df_view import DFView
import polars as pl

class BedpePairs(DFView):
    def __init__(self, df):
        self.df = df

    def reformat(self, fmt):
        from .bedpe_id_pairs import BedpeIDPairs
        BedpePairs.conversions = [BedpeIDPairs]
        self.assert_conversion_allowed(fmt)

        df = self.df.with_row_index("pairID")
        return BedpeIDPairs(df)