import pysam

decompression_threads = 3
compression_threads = 8
write_mode = "w"

def edit_seg1(seg: pysam.AlignedSegment) -> pysam.AlignedSegment:
    RG = seg.get_tag("RG") if seg.has_tag("RG") else ""
    RG = RG[:-2] if RG.endswith(".R") else RG
    seg.set_tag("RG", RG)
    seg.set_tag("ZB", RG)
    return seg

def edit_seg2(seg: pysam.AlignedSegment) -> pysam.AlignedSegment:
    seg.set_tag("ZD", seg.is_duplicate)
    return seg