import duckdb
import polars as pl

def ends_to_intervals(
        pairs, 
        intervals, 
        idx1: str, 
        start1: str, 
        end1: str, 
        idx2: str, 
        start2: str, 
        end2: str, 
        pair_type: bool = True, 
        idx_null =  -1, 
        pos_null = 0,
        memory_limit = None,
        threads = None):
    limits = ""
    if memory_limit: limits += f"SET memory_limit = '{memory_limit}';\n"
    if threads: limits += f"SET threads = {threads};\n"

    sql = f"""
{limits}
WITH intervals1 AS (SELECT * FROM intervals),
     intervals2 AS (SELECT * FROM intervals)

SELECT
    chrom1, 
    pos1,
    chrom2,
    pos2,
    {'pair_type,' if pair_type else ''}
    intervals1.index AS {idx1}, 
    intervals1.start AS {start1},
    intervals1.end AS {end1},
    intervals2.index AS {idx2}, 
    intervals2.start AS {start2},
    intervals2.end AS {end2}

FROM pairs
LEFT JOIN intervals1
ON
    pairs.chrom1 == intervals1.chrom
    AND (pairs.pos1 BETWEEN intervals1.start AND intervals1.end - 1)
LEFT JOIN intervals2
ON pairs.chrom2 == intervals2.chrom
    AND (pairs.pos2 BETWEEN intervals2.start AND intervals2.end - 1)
ORDER BY chrom1, chrom2, pos1, pos2{', pair_type' if pair_type else ''}
"""

    try:
        result = duckdb.execute(
            sql
        ).pl()
    except Exception as e:
        print(f"SQL that generated this error:\n{sql}")
        raise(e)

    try:
        return result.with_columns(
            pl.col(idx1).fill_null(idx_null),
            pl.col(start1).fill_null(pos_null),
            pl.col(end1).fill_null(pos_null),
            pl.col(idx2).fill_null(idx_null),
            pl.col(start2).fill_null(pos_null),
            pl.col(end2).fill_null(pos_null),
        )
    except Exception as e:
        print(f"Failed to fill null columns on ends_to_intervals with idx_null {idx_null} and pos_null {pos_null} after SQL:\n{sql} resulting in dataframe: {result}")
        raise (e)

