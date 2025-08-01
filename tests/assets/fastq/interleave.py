import smart_open
import smart_open_with_pbgzip

lines1 = smart_open.smart_open("CEMBA3C_18B3C_R2_P5-2-M17-L15.R1.fq.gz", "rt").readlines()
lines2 = smart_open.smart_open("CEMBA3C_18B3C_R2_P5-2-M17-L15.R2.fq.gz", "rt").readlines()
count = len(lines1) + len(lines2)
output = open("CEMBA3C_18B3C_R2_P5-2-M17-L15_interleaved.fastq", "w")
l1 = 0
l2 = 0
while l1 < len(lines1) or l2 < len(lines2):
    if l1 <= l2:
        lines = lines1[l1:l1+4]
        l1 += 4
    else:
        lines = lines2[l2:l2+4]
        l2 += 4
    for line in lines:
        output.write(line)
