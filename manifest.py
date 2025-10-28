from hich_manifest.record.records import Record
from hich_manifest.record.entrypoints import Entrypoint
from hich_manifest.config.task_control import TaskControl
from hich_manifest.config.config_fastq_align import ConfigFastqAlign, BWA, BWAMETH, FASTQ_TYPE
from hich_manifest.config.config_bam_parse_pairs import ConfigBamParsePairs
from hich_manifest.config.config_pairs_label import ConfigPairsLabel
from hich_manifest.config.config_pairs_select import ConfigPairsSelect, PairsFilters
from hich_manifest.config.config_pairs_cool_bin import ConfigPairsCoolBin
from hich_manifest.config.command import Command

import pprint
from pydantic import ValidationError

fastq1 = "tests/assets/fastq/1k/1k_ERR1413593_1.fastq.gz"
fastq2 = "tests/assets/fastq/1k/1k_ERR1413593_2.fastq.gz"
aligner_index_dir = "tests/assets/index/bwa"
chromsizes = "tests/assets/chromsizes/M129.sizes"
fragment_index = "tests/assets/fragmentIndex/M129_HindIII.bed"

try:
    record = Record(
        id = "sample1",
        entrypoint = Entrypoint.FASTQ_ALIGN,
        fastq1 = fastq1,
        fastq2 = fastq2,
        aligner_index_dir = aligner_index_dir,
        chromsizes=chromsizes,
        fragment_index = fragment_index,
        bin_resolutions=[1000,2000,5000],
        config_fastq_align = ConfigFastqAlign(fastq_type=FASTQ_TYPE.PAIRED_END, aligner=BWA.MEM2, aligner_index_prefix="M129"),
        config_pairs_select = ConfigPairsSelect(
            pairs_filters=PairsFilters(min_dist_ff=1000, min_dist_fr=1000, keep_pair_chroms=PairsFilters.CisTrans.IS_CIS),
        )
    )
    print(record.model_dump_json(indent=True))

except ValidationError as e:
    pprint.pprint(e.errors())