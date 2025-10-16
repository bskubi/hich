nlines_gz () {
    curl $2 | gunzip -c | head -n $(($1 * 4)) | gzip -c > $3
}

nlines_gz 4000 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR141/003/ERR1413593/ERR1413593_1.fastq.gz 1k_ERR1413593_1.fq.gz
nlines_gz 4000 ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR141/003/ERR1413593/ERR1413593_2.fastq.gz 1k_ERR1413593_2.fq.gz
 
