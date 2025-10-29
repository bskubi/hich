from hich_manifest.record.records import Record, FastqAlignRecord, PairsMergeBeforeDedupTargetRecord, PairsMergeAfterDedupTargetRecord
from hich_manifest.record.entrypoints import Entrypoint
from hich_manifest.config.task_control import TaskControl
from hich_manifest.config.config_fastq_align import ConfigFastqAlign, BWA, BWAMETH, FASTQ_TYPE
from hich_manifest.config.config_bam_parse_pairs import ConfigBamParsePairs
from hich_manifest.config.config_pairs_label import ConfigPairsLabel
from hich_manifest.config.config_pairs_select import ConfigPairsSelect, PairsFilters
from hich_manifest.config.config_pairs_merge import ConfigPairsMerge
from hich_manifest.config.config_pairs_cool_bin import ConfigPairsCoolBin
from hich_manifest.config.command import Command
from hich_manifest.manifest import Manifest
import pprint
from pydantic import ValidationError

fastq1 = "tests/assets/fastq/1k/1k_ERR1413593_1.fastq.gz"
fastq2 = "tests/assets/fastq/1k/1k_ERR1413593_2.fastq.gz"
aligner_index_dir = "tests/assets/index/bwa"
chromsizes = "tests/assets/chromsizes/M129.sizes"
fragment_index = "tests/assets/fragmentIndex/M129_HindIII.bed"

template = FastqAlignRecord.model_construct(**dict(
    entrypoint = Entrypoint.FASTQ_ALIGN,
    aligner_index_dir = aligner_index_dir,
    chromsizes=chromsizes,
    fragment_index = fragment_index,
    bin_resolutions=[1000,2000,5000],
    config_fastq_align = ConfigFastqAlign(fastq_type=FASTQ_TYPE.PAIRED_END, aligner=BWA.MEM2, aligner_index_prefix="M129"),
    config_pairs_select = ConfigPairsSelect(
        pairs_filters=PairsFilters(min_dist_ff=1000, min_dist_fr=1000, keep_pair_chroms=PairsFilters.CisTrans.IS_CIS),
    )
))

try:
    record1 = template.model_copy(update = dict(id="11", fastq1=fastq1, fastq2=fastq2))
    record2 = template.model_copy(update = dict(id="12", fastq1=fastq1, fastq2=fastq2))
    record3 = Record(
        id = "1_before", 
        chromsizes=chromsizes, 
        entrypoint=Entrypoint.PAIRS_MERGE_BEFORE_DEDUP_TARGET, 
        config_pairs_merge = ConfigPairsMerge(source_ids=["11", "12"])
    )
    record4 = Record(
        id = "1_after",
        chromsizes=chromsizes,
        entrypoint=Entrypoint.PAIRS_MERGE_AFTER_DEDUP_TARGET,
        config_pairs_merge = ConfigPairsMerge(source_ids=["11", "12", "1_after"])
    )
    
    manifest = Manifest()
    manifest.set_template(template)
    manifest.add_record(dict(id="11", fastq1=fastq1, fastq2=fastq2))
    manifest.add_record(dict(id="12", fastq1=fastq1, fastq2=fastq2))
    manifest.add_record(
        PairsMergeBeforeDedupTargetRecord(
            id = "1_before", 
            chromsizes=chromsizes, 
            entrypoint=Entrypoint.PAIRS_MERGE_BEFORE_DEDUP_TARGET, 
            config_pairs_merge = ConfigPairsMerge(source_ids=["11", "12"])
        )
    )
    manifest.add_record(
        PairsMergeAfterDedupTargetRecord(
            id = "1_after",
            chromsizes=chromsizes,
            entrypoint=Entrypoint.PAIRS_MERGE_AFTER_DEDUP_TARGET,
            config_pairs_merge = ConfigPairsMerge(source_ids=["11", "12", "1_after"])
        )
    )
    with open("manifest.json", "w") as file:
        file.write(manifest.model_dump_json(indent=True))
    
except ValidationError as e:
    pprint.pprint(e.errors())