from hich_manifest.fastq_align import FastqAlign, BWA
import pprint
from pydantic import ValidationError

fastq1 = "tests/assets/fastq/1k/1k_ERR1413593_1.fastq.gz"
fastq2 = "tests/assets/fastq/1k/1k_ERR1413593_2.fastq.gz"
aligner_index_dir = "tests/assets/index/bwa"
chromsizes = "tests/assets/chromsizes/M129.sizes"
fragment_index = "tests/assets/fragmentIndex/M129_HindIII.bed"
try:
    print(FastqAlign(fastq1=fastq1, fastq2=fastq2, aligner=BWA.MEM2, aligner_index_dir=aligner_index_dir, aligner_index_prefix="M129"))
except ValidationError as e:
    pprint.pprint(e.errors())