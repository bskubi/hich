import pysam

decompression_threads = 3
compression_threads = 8
write_mode = "w"

def edit_seg1(seg: pysam.AlignedSegment) -> pysam.AlignedSegment:
    return seg

def edit_seg2(seg: pysam.AlignedSegment) -> pysam.AlignedSegment:
    seg.set_tag("ZD", seg.is_duplicate)
    return seg