from enum import Enum
from pathlib import Path
from typing import Optional

from pydantic import model_validator, Field, BaseModel, FilePath
from .config_process_task_control import ConfigProcessTaskControl
from .command import Command



class PairsFilters(BaseModel):
    class CisTrans(str, Enum):
        IS_CIS = "CIS"
        IS_TRANS = "TRANS"
        BOTH = "BOTH"
    keep_pair_types: Optional[list[str]] = Field(["UU", "UR", "RU"], descripton="Pair end types to keep (U=unique alignment, R=rescue of unobserved ligation junction)")
    min_dist_ff: Optional[int] = Field(None, ge=1, description="Minimum cis end-to-end distance for ++ pairs")
    min_dist_fr: Optional[int] = Field(None, ge=1, description="Minimum cis end-to-end distance for +- pairs")
    min_dist_rf: Optional[int] = Field(None, ge=1, description="Minimum cis end-to-end distance for -+ pairs")
    min_dist_rr: Optional[int] = Field(None, ge=1, description="Minimum cis end-to-end distance for -- pairs")
    keep_pair_chroms: CisTrans = Field(CisTrans.BOTH, description="Keep cis (intrachromosomal) pairs, trans (interchromosomal) pairs, or both")
    not_same_frag: bool = Field(True, description="Drop pairs that map to same restriction fragment")

    @property
    def CONDITION(self) -> str:
        conditions = []
        
        # Pair alignment types (i.e. U = unique, R = rescue)
        if self.keep_pair_types:
            conditions.append(f"pair_type in {self.keep_pair_types}")
        
        # Cis/trans filters
        if self.keep_pair_chroms == PairsFilters.CisTrans.IS_TRANS:
            conditions.append("chrom1 != chrom2")
        elif self.keep_pair_chroms == PairsFilters.CisTrans.IS_CIS:
            conditions.append("chrom1 == chrom2")
        
        # Orientation and distance filters
        if any([self.min_dist_ff, self.min_dist_fr, self.min_dist_rf, self.min_dist_rr]):
            orient_dist_filters = []
            if self.min_dist_ff:
                orient_dist_filters.append(f"(strand1+strand2=='++' and pos2-pos1 >= {self.min_dist_ff})")
            if self.min_dist_fr:
                orient_dist_filters.append(f"(strand1+strand2=='+-' and pos2-pos1 >= {self.min_dist_fr})")
            if self.min_dist_rf:
                orient_dist_filters.append(f"(strand1+strand2=='--' and pos2-pos1 >= {self.min_dist_rf})")
            if self.min_dist_rr:
                orient_dist_filters.append(f"(strand1+strand2=='--' and pos2-pos1 >= {self.min_dist_rr})")
            if orient_dist_filters:
                orient_dist_filters = ["chrom1 != chrom2"] + orient_dist_filters
                conditions.append(" or ".join(orient_dist_filters))

        # Ends map to different restriction fragments
        if self.not_same_frag:
            conditions.append("rfrag1 != rfrag2")

        # Require all conditions to pass
        CONDITION = " and ".join([f"({c})" for c in conditions])
        return CONDITION


class ConfigPairsSelect(ConfigProcessTaskControl):
    pairs_filters: PairsFilters = PairsFilters()
    chrom_subset: Optional[FilePath] = None

    command_pairtools_select: Optional[Command] = Field(
        default_factory = lambda data: Command(
                base_command = "pairtools select",
                opts = {
                    "--output": "${pairs_selected}",
                    "--nproc-in": "${cpus}",
                    "--nproc-out": "${cpus}",
                    **({"--chrom-subset": data.get('chrom_subset').name} if data.get('chrom_subset') is not None else {})
                },
                args = [
                    data['pairs_filters'].CONDITION if data.get('pairs_filters') else "true", 
                    "${pairs}"
                ]
            )
    )
    