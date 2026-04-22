TODO:
0. Add consolidation and vacuum to Hich engine
1. Remove limit-chroms in write DNA methylation section.
2. Validate download.sh and index.sh test setup scripts
3. Allow filtering output write and statistics input separately
4. Test: interleaved, single-end alignment
5. Test: aligner combinations
6. Restriction fragment filtering
7. Dockerize, distribute
8. Test suite

```
cat input.bed | \
sort -k1,1 -k2,2n | \
bedtools genomecov -bg -i stdin -g chrom.sizes | \
bedGraphToBigWig stdin chrom.sizes output.bw
```