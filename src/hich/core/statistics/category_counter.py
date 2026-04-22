import polars as pl
from dataclasses import dataclass, field
from collections import Counter
from typing import Callable, Optional, List, Generic, TypeVar
import ast
from .statistics import Statistic, StatisticType, StatisticAccumulator
from loguru import logger

@dataclass
class Category(StatisticType):
    """

    Attributes:
        stats (Optional[Callable[[pl.DataFrame], pl.DataFrame]]): Function taking DataFrame raw observations and returning dataframe of category observations 
        pivot_on (Optional[List[str]]): If given, final output will utilize these category dimensions as pivot columns.
        pivot_index (Optional[List[str]]): If given, final outptu will utilize these category dimensions as pivot index.
        count_col (str): Temporary name of category count column.
        category_dims: Return the columns of 'df'
        category_counts: Return the number of observations for each category in self.df

    Example:
        ```python
        >>> import polars as pl
        >>> from hich.engine.category_count import CategoryCount
        >>> def get_stats(df: pl.DataFrame) -> pl.DataFrame:
                return df.select("colA")
        >>> category_count = CategoryCount(get_stats)
        >>> obs_df = pl.DataFrame({"colA":[1,1,2], "colB":[1,2,3]})
        >>> category_count.stats(obs_df)
        shape: (3, 1)
        ┌──────┐
        │ colA │
        │ ---  │
        │ i64  │
        ╞══════╡
        │ 1    │
        │ 1    │
        │ 2    │
        └──────┘
        ```
    """
    categorize: Optional[Callable[[pl.DataFrame], pl.DataFrame]] = None
    pivot_on: Optional[List[str]] = None
    pivot_index: Optional[List[str]] = None
    rename_pivot_on_pattern: Optional[str] = None
    count_column: str = "CategoryCount"

    @property
    def accumulator(self) -> type["CategoryCounter"]:
        return CategoryCounter

@dataclass
class CategoryCounter(StatisticAccumulator):
    """
    Manage batched extraction of category counts.

    Attributes:
        category (Category): Defines how to extract categories from observations and format category counts DataFrame
        counter (Counter): Current counts for each category
        schema (Optional[pl.Schema]): Schema of categories DataFrame 
    """
    category: Category
    counter: Counter = field(default_factory=Counter)
    schema: Optional[pl.Schema] = None

    @property
    def statistic(self) -> Statistic:
        """
        Wrap 'self.counts' as a Statistic object

        Returns:
            Statistic: self.counts
        """
        return Statistic(self.counts)

    def update(self, obs_df: pl.DataFrame):
        """
        Update category counts with new observations.

        Arguments:
            obs_df (pl.DataFrame): Observations to extract categories from.
        """
        categories_df = self._categories(obs_df)
        self.set_schema(categories_df)
        self.counter.update(df_to_counter(categories_df))

    def set_schema(self, categories_df: pl.DataFrame):
        """
        Store schema of categories_df.

        Raises:
            ValueError: schema is set but does not match schema of categories_df
        """
        # Validate category dim consistency
        if self.schema is not None:
            category_dims_match = self.schema == categories_df.schema
            if not category_dims_match:
                err_msg = f"Schema mismatch:\n{categories_df.schema}\n{categories_df}\n{self}"
                raise ValueError(err_msg)

        self.schema = categories_df.schema

    @property
    def counter_to_df(self) -> pl.DataFrame:
        return counter_to_df(self.counter, self.schema, self.category.count_column)

    @property
    def counts(self) -> pl.DataFrame:
        """
        Convert counter and schema to Polars DataFrame

        Attributes:
            counter (Counter): The counter to convert
            schema (pl.Schema): The schema of the output DataFrame not including the count column
        """
        # Don't pivot, just return as-is
        if not self.category.pivot_on:
            return self.counter_to_df
        
        return pivot_and_rename(
            df=self.counter_to_df, 
            on=self.category.pivot_on, 
            index=self.category.pivot_index, 
            values=self.category.count_column, 
            pattern=self.category.rename_pivot_on_pattern
        )
        

    def _categories(self, obs_df: Optional[pl.DataFrame] = None) -> Optional[pl.DataFrame]:
        """
        Extract categories from observations by calling self.categorize(obs_df)

        Arguments:
            obs_df (Optional[pl.DataFrame]): Observations DataFrame passed to self.categorize().

        Returns:
            pl.DataFrame: categories DataFrame (or None if obs_df is None)
        """
        return obs_df if self.category.categorize is None or obs_df is None else self.category.categorize(obs_df)    

@logger.catch(reraise=True)
def df_to_counter(df: Optional[pl.DataFrame]) -> Counter:
    """
    Count unique rows in DataFrame

    Attributes:
        df (pl.DataFrame): Number of observations for each unique row are summed
        count_column (str): Temporary count column to add to dataframe
    
    self.schema
    """
    if df is None:
        return Counter()
    
    # Create guaranteed unique temporary count column
    count_column = "__" + "".join(df.columns) + "__"

    # Count observations of categories
    stats_df = (
        df
        .group_by(*df.columns)
        .agg(pl.len().alias(count_column))
    )
    
    # Convert categories to list of tuples
    categories = list(stats_df.select(*df.columns).iter_rows())
    
    # Get counts for each category
    counts = stats_df[count_column]
    
    # Create Counter object
    category_counts_dict = {category: count for category, count in zip(categories, counts)}
    counter = Counter(category_counts_dict)
    return counter

@logger.catch(reraise=True)
def counter_to_df(counter: Counter, schema: pl.Schema, count_column: str) -> pl.DataFrame:
    """
    Convert Counter to polars DataFrame, putting count values into column with alias count_column

    Arguments:
        counter (Counter): A Counter argument containing columns matching schema
        schema (pl.Schema): A Polars Schema matching the counter keys
        count_column (str): Alias for the count values in the counter

    Returns:
        pl.DataFrame: Rows match keys and count in the counter, schema is 'schema' with extra UInt64 column for count_column. 
    """
    if not counter:
        count = pl.Series([], dtype=pl.UInt64)
        
        df = pl.DataFrame(schema=schema)
    else:
        count = pl.Series(counter.values(), dtype=pl.UInt64)
        df = pl.from_records(list(counter.keys()), schema=schema, orient='row')
    
    df = df.with_columns(count.alias(count_column)).sort(schema.names())
    
    return df

@logger.catch(reraise=True)
def pivot_and_rename(
        df: pl.DataFrame,
        on: str | List[str], 
        index: str | List[str] | None = None,
        values: str | List[str] | None = None, 
        pattern:  Optional[str] = None
    ) -> pl.DataFrame:
    """
    Pivot df and optionally rename pivot columns

    Arguments:
        df (pl.DataFrame): DataFrame to pivot
        on (List[str]): List of pivot column names

    """
    # For pivot:
    # If 'on' + 'values' is proper subset of df.columns, df.columns - on - values is used as index
    # If 'on' + 'index' is proper subset of df.columns, df.columns - on - index is used as values
    # If 'on' + 'values' equals df.columns, raises error
    # 'values' or 'index' must be specified
    # We introduce one extra piece of magic: if 'index' unspecified and 'on' + 'values' equals df.columns,
    # then we produce a single-row dataframe with one entry for each combination of 'on' values.

    # Check if 'on' + 'values' account for all columns.
    columns_set = set(df.columns)
    match on:
        case str(): on_set = set([on])
        case None: on_set = set()
        case _:
            assert all(isinstance(col, str) for col in on), "on must contain exclusively objects of type str"
            on_set = set(on)
    match values:
        case str(): values_set = set([values])
        case None: values_set = set()
        case _: 
            assert all(isinstance(col, str) for col in index), "values must contain exclusively objects of type str"
            values_set = set(values)
    
    # Contains all columns in df.columns that aren't going to be found by polars as in 'on' or 'values'
    index_set = columns_set - on_set - values_set

    if index is None and not index_set:
        # Produce single-row DataFrame instead of raising an exception
        # Example:
        """
        Original df:
        A Count
        X 1
        Y 2

        Add column:
        __ACount__ A Count
        0          X 1
        1          Y 2

        Pivot:
        __ACount__ X    Y
        0          1    null
        1          null 2

        Drop index:
        X    Y
        1    null
        null 2

        Fill null:
        X   Y
        1   2
        1   2

        Unique:
        X   Y
        1   2
        """
        # Create temporary column (i.e. __ACount__)
        temp_col = "__" + "".join(df.columns) + "__"

        # Pivot, fill nulls with max (which just uses the unique non-null value in the column)
        # then take 'unique' (which reduces to a single row)
        pivot = (
            df.with_row_index(temp_col)
            .pivot(on=on, values=values, aggregate_function="sum")
            .drop(temp_col)    
            .fill_null(strategy="max")
            .unique()
        )
    else:
        # Just do a normal pivot
        pivot = df.pivot(on=on, index=index, values=values)

    # Get all columns produced by pivot
    pivot_on_cols = [col for col in pivot.columns if col not in index_set]

    # Replace null columns with 0
    pivot = pivot.with_columns(*[
        pl.col(col).fill_null(0)
        for col in pivot_on_cols
    ])

    if pattern:
        # Rename 'on' columns
        return rename_pivot_on_columns(pivot, on, pattern, index_set)
    else:
        return pivot

@logger.catch(reraise=True)
def rename_pivot_on_columns(pivot_df: pl.DataFrame, on: list[str], pattern: str, ignore: Optional[List[str]]) -> pl.DataFrame:
    """
    Improve formatting of auto-assigned polars pivot column names when multiple columns are used for 'on' argument of pivot method.

    Arguments:
        pivot_df (pl.DataFrame): Polars DataFrame from call to pivot method
        on (list[str]): Columns used in 'on' argument of pivot method.
        pattern (str): f-string pattern used to rename pivot columns. Use column names in 'on' as keys.
        ignore (list[str]): Column names to NOT rename in the pivot dataframe (i.e. columns in 'index' argument to pivot).

    Returns:
        pl.DataFrame: pivot table with renamed columns

    When multiple col names are passed to 'on' argument of Polars DataFrame.pivot,
    the column name generate is formatted with pairs of values, like {"colA_val1, colB_val1"}.
    Multiple index columns are retained as the original columns.

    >>> df = pl.DataFrame(dict(A=[1,2,3], B=[4,5,6], C=[7,8,9], D=[10,11,12], E=[13,14,15]))
    >>> df.pivot(on=["A","B"], index=["C","D"], values="E")
    shape: (3, 5)
    ┌─────┬─────┬───────┬───────┬───────┐
    │ C   ┆ D   ┆ {1,4} ┆ {2,5} ┆ {3,6} │
    │ --- ┆ --- ┆ ---   ┆ ---   ┆ ---   │
    │ i64 ┆ i64 ┆ i64   ┆ i64   ┆ i64   │
    ╞═════╪═════╪═══════╪═══════╪═══════╡
    │ 7   ┆ 10  ┆ 13    ┆ null  ┆ null  │
    │ 8   ┆ 11  ┆ null  ┆ 14    ┆ null  │
    │ 9   ┆ 12  ┆ null  ┆ null  ┆ 15    │
    └─────┴─────┴───────┴───────┴───────┘

    >>> format_pivot_column_names(df, on=["A", "B"], pattern=["{A}_{B}"], ignore=["C","D"])
    ┌─────┬─────┬───────┬───────┬───────┐
    │ C   ┆ D   ┆  1_4  ┆  2_5  ┆  3_6  │
    │ --- ┆ --- ┆ ---   ┆ ---   ┆ ---   │
    │ i64 ┆ i64 ┆ i64   ┆ i64   ┆ i64   │
    ╞═════╪═════╪═══════╪═══════╪═══════╡
    │ 7   ┆ 10  ┆ 13    ┆ null  ┆ null  │
    │ 8   ┆ 11  ┆ null  ┆ 14    ┆ null  │
    │ 9   ┆ 12  ┆ null  ┆ null  ┆ 15    │
    └─────┴─────┴───────┴───────┴───────┘
    """
    print(on, ignore, pivot_df.columns)
    ignore = ignore or []
    use_columns = set(pivot_df.columns).difference(ignore)
    rename = {}
    for col in use_columns:
        if len(on) > 1:
            # Convert from set to list to preserve order on eval
            # Example: '{"1", "2"}' -> '["1", "2"]'
            elements = "[" + col.strip("{}") + "]"

            # Evaluate string representation of list to actual list of column names
            # Converts '["1", "2"]' (str) to ["1", "2"] (list)
            elements = ast.literal_eval(elements)
        else:
            # Single-column pivot, so just store actual column name
            elements = [col]
        
        # Reshape to dict with non-pivot 'on' column names as keys, corresponding values in pivot column as values 
        # ["1", "2"] -> {"A": "1", "B": "2"}
        format_kwargs = {on_col: element for on_col, element in zip(on, elements)}

        # Add to renaming dict. If pattern is '{A}_{B}':
        # rename['{"1", "2"}'] = '{A}_{B}'.format(A = "1", B = "2") = '1_2'
        rename[col] = pattern.format(**format_kwargs)
    
    # Update the pivot column names
    return pivot_df.rename(rename)