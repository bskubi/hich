import duckdb
import polars as pl
import numpy as np





def map_ends(
        pairs: pl.DataFrame, 
        intervals: pl.DataFrame, 
        pairs_idx1: str,
        pairs_pos1: str,
        pairs_idx2: str,
        pairs_pos2: str,
        ivl_idx:str,
        ivl_start: str,
        ivl_end: str,
        tag_idx1: str, 
        tag_start1: str, 
        tag_end1: str, 
        tag_idx2: str, 
        tag_start2: str, 
        tag_end2: str,
        idx_null =  -1, 
        pos_null = 0,
        ridx = "__row_index__"):
    
    columns = pairs.columns + [tag_idx1, tag_start1, tag_end1, tag_idx2, tag_start2, tag_end2]

    assert not ridx in pairs, f"Pairs file must not contain column named {ridx}, but received chunk:\n{pairs}"
    assert not ridx in intervals, f"Intervals file must not contain column named {ridx}, but received:\n{intervals}"
    pairs = pairs.with_row_index(ridx)

    def unique_pairs_vals(col):
        vals = pairs[col]
        vals = vals.unique()
        vals = vals.to_list()
        return vals

    idx1_unique = unique_pairs_vals(pairs_idx1)
    idx2_unique = unique_pairs_vals(pairs_idx2)
    indices = idx1_unique + idx2_unique
    indices = sorted(list(set(indices)))
    
    
    tag1_chunks = []
    tag2_chunks = []

    for idx in indices:
        try:
            chunk1 = pairs.drop(pairs_idx2, pairs_pos2).filter(pl.col(pairs_idx1) == idx)
            chunk1 = chunk1.sort(pairs_pos1)
            chunk2 = pairs.select(ridx, pairs_idx2, pairs_pos2).filter(pl.col(pairs_idx2) == idx)
            chunk2 = chunk2.sort(pairs_pos2)
            chunk_intervals = intervals.filter(pl.col(ivl_idx) == idx)
            chunk_intervals = chunk_intervals.with_row_index(ridx)
        except Exception as e:
            print(f"Failed to sort on {pairs_idx1} and {pairs_idx2}. Pairs:\n{pairs}")
            raise(e)

        if idx not in chunk_intervals[ivl_idx]:
            tag1_chunk = chunk1.with_columns(
                pl.lit(idx_null).alias(tag_idx1).cast(pl.String),
                pl.lit(pos_null).alias(tag_start1),
                pl.lit(pos_null).alias(tag_end1)
            )
            tag1_chunks.append(tag1_chunk)

            tag2_chunk = pl.DataFrame({
                ridx: chunk2[ridx],
                pairs_idx2: chunk2[pairs_idx2],
                pairs_pos2: chunk2[pairs_pos2],
                tag_idx2: pl.Series([idx_null]*len(chunk2)).cast(pl.String),
                tag_start2: [pos_null]*len(chunk2),
                tag_end2: [pos_null]*len(chunk2)
            })
            tag2_chunks.append(tag2_chunk)
        else:
            idx_intervals = chunk_intervals.filter(pl.col(ivl_idx) == idx)
            try:
                indices = idx_intervals[ridx]
                starts = idx_intervals[ivl_start]
                ends = idx_intervals[ivl_end]
            except Exception as e:
                print(f"Failed to find one of {[ivl_idx, ivl_start, ivl_end]} in columns of intervals:\n{idx_intervals}")

            chunk1_ssidx = np.searchsorted(ends, chunk1[pairs_pos1], 'right')
            chunk2_ssidx = np.searchsorted(ends, chunk2[pairs_pos2], 'right')

            assert len(chunk1[pairs_pos1]) == 0 or max(chunk1_ssidx) < len(indices), f"{pairs_pos1} had value {max(chunk1[pairs_pos1])} >= maximum end position {max(ends)}"
            assert len(chunk2[pairs_pos2]) == 0 or max(chunk2_ssidx) < len(indices), f"{pairs_pos2} had value {max(chunk2[pairs_pos2])} >= maximum end position {max(ends)}"
            
            tag1_idx = indices[chunk1_ssidx]
            tag1_start = starts[chunk1_ssidx]
            tag1_end = ends[chunk1_ssidx]

            tag1_chunk = chunk1.with_columns(
                pl.Series(tag1_idx).alias(tag_idx1).cast(pl.String),
                pl.Series(tag1_start).alias(tag_start1),
                pl.Series(tag1_end).alias(tag_end1),
            )
            tag1_chunks.append(tag1_chunk)

            tag2_ridx = chunk2[ridx]
            tag2_idx = indices[chunk2_ssidx]
            tag2_start = starts[chunk2_ssidx]
            tag2_end = ends[chunk2_ssidx]
            tag2_chunk = pl.DataFrame({
                ridx: tag2_ridx,
                pairs_idx2: chunk2[pairs_idx2],
                pairs_pos2: chunk2[pairs_pos2],
                tag_idx2: pl.Series(tag2_idx).cast(pl.String),
                tag_start2: tag2_start,
                tag_end2: tag2_end
            })
            tag2_chunks.append(tag2_chunk)

    tag1_chunks = pl.concat(tag1_chunks, how = 'vertical_relaxed')
    tag1_chunks = tag1_chunks.sort(ridx)
    tag2_chunks = pl.concat(tag2_chunks, how = 'vertical_relaxed')
    tag2_chunks = tag2_chunks.sort(ridx)
    tag1_chunks = tag1_chunks.drop(ridx)
    tag2_chunks = tag2_chunks.drop(ridx)
    tagged_chunks = pl.concat([tag1_chunks, tag2_chunks], how = 'horizontal')
    tagged_chunks = tagged_chunks.select(columns)
    return tagged_chunks
