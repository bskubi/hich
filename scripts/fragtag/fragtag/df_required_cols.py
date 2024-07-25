from .df_view import DFView
import polars as pl

class DFRequiredCols(DFView):
    def __init__(self, df):
        super().__init__(df)
        
        # Ensure all required columns are present in the dataframe
        required_cols = set(type(self).required_cols)
        cols = set(self.df.columns)
        missing = required_cols.difference(cols)

        err = f"{type(self)} instance is missing columns {missing}"
        assert not missing, err
        