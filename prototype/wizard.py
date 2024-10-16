from multipledispatch import dispatch
import param
import copy
from statemachine import StateMachine, State
import subprocess
from yaml import CLoader as Loader, CDumper as Dumper
import yaml

def load_nextflow_config(filename):
    yamlstr = subprocess.run(f"groovy toyaml.groovy {filename}", shell=True, stdout = subprocess.PIPE).stdout
    object = yaml.load(yamlstr, Loader = Loader)
    return object

class _Blank: pass
Blank = _Blank()

class _Comment:
    def __init__(self, min_split = 80, level_size = 4):
        self.min_split = min_split
        self.level_size = 4
    
    def format(self, comment_str, level):
        lines = []
        tabs = "\t"*level
        comment_str = comment_str.split()
        current_line = tabs + "//"
        for word in comment_str:
            current_line += " " + word
            row_len = level*(self.level_size - 1) + len(current_line) + len(word) + 1
            if row_len > self.min_split:
                lines.append(current_line)
                current_line = tabs + "//"
        if current_line != tabs + "//":
            lines.append(current_line)
        result = "\n".join(lines)
        return result

class BracedDict(dict): pass

Comment = _Comment()

@dispatch(str, bool, int)
def to_nextflow_config(key, bool_value, level):
    tabs = "\t"*level
    return f"{tabs}{key} = " + ("true" if bool_value else "false")

@dispatch(str, int)
def to_nextflow_config(string_value, level):
    tabs = "\t"*level
    return tabs + string_value

@dispatch(str)
def to_nextflow_config(str_value):
    return f'"{str_value}"'

@dispatch(object)
def to_nextflow_config(object_value):
    return "null" if object_value is None else str(object_value)

@dispatch(str, _Blank, int)
def to_nextflow_config(key, blank, level):
    tabs = "\t"*level
    return f"{tabs}{key}"

@dispatch(str, _Comment, int)
def to_nextflow_config(key, comment_value, level):
    return comment_value.format(key, level)

@dispatch(str, str, int)
def to_nextflow_config(key, str_value, level):
    tabs = "\t"*level
    
    return f'{tabs}{key} = "{str_value}"'

@dispatch(str, object, int)
def to_nextflow_config(key, object_value, level):
    tabs = "\t"*level
    if object_value is None:
        object_value = "null"
    return f"{tabs}{key} = {object_value}"

@dispatch(str, dict, int)
def to_nextflow_config(key, dict_value, level):
    tabs = "\t"*level
    if any([isinstance(value, dict) for value in dict_value.values()]) or len(dict_value) > 3 or isinstance(dict_value, BracedDict):
        start = f"{tabs}{key} {{\n"
        middle = []
        for key, value in dict_value.items():
            middle.append(to_nextflow_config(key, value, level + 1))
        end = f"\n{tabs}}}"

        result = start + "\n".join(middle) + end
        return result
    else:
        key_value_pairs = [": ".join([str(key), to_nextflow_config(value)]) for key, value in dict_value.items()]
        hashmap_entries = ", ".join(key_value_pairs)
        hashmap = f"[{hashmap_entries}]"
        result = f"{tabs}{key} = {hashmap}"
        return result

@dispatch(dict)
def to_nextflow_config(dict_value):
    result = []
    for key, value in dict_value.items():
        result.append(to_nextflow_config(key, value, 0))
    return "\n".join(result)

config = {
    "singularity.enabled": True,
    "params": {
        "general": {
            "sampleFileSep": "\\t",
            "humid": None,
            "publish": {
                "Nextflow publishDir param for all processes": Comment,
                "See: https://www.nextflow.io/docs/latest/process.html#publishdir": Comment,
                "mode": "copy",
                "genomeReference": "resources/hich/genomeReference",
                "chromsizes": "resources/hich/chromsizes",
                "bwaMem2Index": "resources/hich/bwa-mem2",
                "bwaMemIndex": "resources/hich/bwa-mem",
                "fragmentIndex": "resources/hich/fragmentIndex",
                "align": "results/align",
                "parse": "results/pairs/parse",
                "dedup": "results/pairs/dedup",
                "mcool": "results/matrix/mcool",
                "hic": "results/matrix/hic",
                "pairStats": "results/pairStats",
                "qc": "results/qc"
            },
            """After these steps, generate read-level QC stats and present in a combined
            MultiQC report containing all samples at each processing stage.""": Comment,
            "qcAfter": ["parse", "ingestPairs", "tagRestrictionFragments", "deduplicate", "select"],
            """Number of reads to downsample to (by default) for a humid run, which can also be modified with --humid [numberOfReads]""": Comment,
            """Example: '--humid 400000'""": Comment,
            "humidDefault": 100000
        },

        "sampleSelectionStrategies": {
            "fairComparisons": BracedDict({
                "same": ["aggregateProfileName", "aggregateLevel"],
                "aggregateProfileName": "coverageMatched"
            }),
            "conditionsOnly": {
                "aggregateLevel": "condition"
            },
            "fullCoverage": {
                "aggregateProfileName": "fullCoverage"
            }
        },

        "aggregate": BracedDict({}),

        "defaults": {
            "techrep": 1,
            "biorep": 1,
            "minMapq": 1,

            "aligner": "bwa-mem2",
            "bwaFlags": "-SP5M",
            "reshapeParams": [],
            "globalDefaultReshapetoCellID": {"option":"--regex", "pattern":"^(.*?):.*", "group":1},

            "pairsFormat": {"chrom1": 2, "pos1": 3, "chrom2": 4, "pos2": 5},

            "parseParams": ["--flip", "--min-mapq 1", "--drop-readid", "--drop-seq", "--drop-sam"],

            "selectFilters": {
                "keepPairTypes": ["UU", "RU", "UR"],
                "keepTrans": True,
                "keepCis": True,
                "minDistFR": 1000,
                "minDistRF": 1000,
                "minDistFF": 1000,
                "minDistRR": 1000,
                "discardSingleFrag": True
            },

            "coolerZoomifyParams": ["--balance", "--balance-args '--max-iters 2000 --trans-only'"],

            "matrix": BracedDict({
                "makeMcoolFileFormat": True,
                "makeHicFileFormat": True,
                "resolutions": [25000, 500000, 1000000, 3000000]
            })
        }
    }
}

default_profiles = {
    "localPC": {
        "executor.name": "local",
        "executor.cpus": 10,
        "executor.memory": "20.GB",
        "process.executor": "local",

        "process": {
            "withLabel: whenLocal_allConsuming": {
                "maxForks": 1,
                "cpus": 10,
                "memory": "20.GB"
            },

            "withLabel: convertHicToMcool": {
                "cpus": 2
            },

            "withLabel: convertMcoolToHic": {
                "cpus": 10
            }
        }
    },

    "localHPC": {
        "executor.name": "local",
        "executor.cpus": 100,
        "executor.memory": "200.GB",
        "process.executor": "local",
        "process": {
            "withLabel: whenLocal_allConsuming": {
                "maxForks": 1,
                "cpus": 100,
                "memory": "200.GB"
            },

            "withLabel: convertMcoolToHic": {
                "cpus": 20
            }
        }
    },

    "jobArray": {
        "process": {
            "withLabel: doJobArray": {
                "array": 20
            }
        }
    },

    "grid": {
        "process": {
            "withLabel: align": {
                "time": "36hr",
                "cpus": 24
            },

            "withLabel: index": {
                "time": "36hr",
                "cpus": 24
            },

            "withLabel: smallResource": {
                "time": "4hr",
                "cpus": 8,
                "memory": "8.GB"
            },

            "withLabel: pairs": {
                "time": "16hr",
                "cpus": 8,
                "memory": "8.GB"
            },

            "withLabel: createMatrix": {
                "time": "36hr",
                "cpus": 8,
                "memory": "16.GB"
            },

            "withLabel: convertHicToMcool": {
                "cpus": 2,
                "memory": "8.GB"
            },

            "withLabel: convertMcoolToHic": {
                "cpus": 20,
                "memory": "8.GB"
            }
        }
    }
}

class ParameterizedConfig: pass

@dispatch(str, ParameterizedConfig, int)
def to_nextflow_config(key, config_value, level):
    return to_nextflow_config(key, config_value.to_dict(), level)

@dispatch(ParameterizedConfig)
def to_nextflow_config(config_value):
    return to_nextflow_config("", config_value.to_dict(), 0)

class AggregateProfile(param.Parameterized, ParameterizedConfig):
    dedupSingleCell = param.Boolean(
        False,
        doc = "Only treat reads as duplicates if they originate from the same cell (have matching cellID column values)"
    )
    dedupMethod = param.Selector(
        objects = ["max", "sum"],
        instantiate = True,
        doc = (
        "During deduplication, reads are ordered by (chrom1, chrom2, pos1, pos2). "
        "Two consecutive reads (readA and readB) are compared. A distance metric "
        "for readA and readB is computed determined by dedupMethod being 'max' or 'sum'. "
        "A necessary but not always sufficient condition for readB to be called as a duplicate "
        "of readA is that this distance metric is at or below a threshold specified by dedupMaxMismatch. "
        "'max' means the distance metric is the maximum distance between readA.pos1 and readB.pos1 and readA.pos2 and readB.pos2. "
        "'sum' means that the sum of these two distances is used."
        )
    )
    dedupMaxMismatch = param.Integer(
        3,
        doc = (
            "For treating reads as duplicates, use this as the maximum threshold for the distance "
            "metric computed via dedupMethod for two consecutive reads to be potentially called as duplicates."
        )
    )
    techrepDedup = param.Boolean(
        True,
        doc = "After the (optional) techrep->biorep merge, deduplicate the techrep-level samples"
    )
    biorepDedup = param.Boolean(
        True,
        doc = "Prior to the (optional) biorep->condition merge, deduplicate the biorep-level samples"
    )
    conditionDedup = param.Boolean(
        False,
        doc = "Deduplicate the condition-level samples (defaults to False to avoid treating linearly proximal reads from separate cells as PCR duplicates)"
    )
    mergeTechrepToBiorep = param.Boolean(True, doc = "Merge techrep-level samples from the same condition and biorep into biorep-level samples")
    mergeBiorepToCondition = param.Boolean(True, doc = "Merge biorep-level samples from the same condition into condition-level samples")
    techrepReadConjuncts = param.ListSelector(objects = ["chrom1", "chrom2", "abs(pos1-pos2)", "cellID", "(rfrag1 == rfrag2)", "pair_type", "stratum"],
    check_on_set = False,
    doc = ("For statistics calling and coverage homogenization, Hich partitions the reads from each sample into blocks defined by the intersection of a set of "
           "CONJUNCTS, which are traits of the read. Reads having the same value for all CONJUNCTS will be assigned to the same READ BLOCK, and you can "
           "have Hich group your samples and downsample them to have the largest possible coverage while matching coverage across all READ BLOCKS. "
           "The benefit of this process is that it allows you to control for coverage differences across samples, as might result from different sequencing depths "
           "or techical errors during the assay. The techrepReadConjuncts parameter sets which CONJUNCTS should be used for this downsampling-based coverage control process "
           "specifically for techrep-level samples."))
    techrepCisStrata = param.List(
        [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000],
        item_type = int,
        doc = ("Partition the cis (intrachromosomal) TECHREP-level contacts into strata with these insert size boundaries "
               "and homogenize coverage within each stratum. Note that these are not the binning resolutions for the contact matrix. "
               "They are potential conjuncts used to homogenize the distribution of contacts during aggregation, prior to binning as a contact matrix.")
    )
    biorepReadConjuncts = param.ListSelector(objects = ["chrom1", "chrom2", "abs(pos1-pos2)", "cellID", "(rfrag1 == rfrag2)", "pair_type", "stratum"],
    check_on_set = False,
    doc = ("For statistics calling and coverage homogenization, Hich partitions the reads from each sample into blocks defined by the intersection of a set of "
           "CONJUNCTS, which are traits of the read. Reads having the same value for all CONJUNCTS will be assigned to the same READ BLOCK, and you can "
           "have Hich group your samples and downsample them to have the largest possible coverage while matching coverage across all READ BLOCKS. "
           "The benefit of this process is that it allows you to control for coverage differences across samples, as might result from different sequencing depths "
           "or techical errors during the assay. The biorepReadConjuncts parameter sets which CONJUNCTS should be used for this downsampling-based coverage control process "
           "specifically for biorep-level samples."))
    biorepCisStrata = param.List(
        [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000],
        item_type = int,
        doc = ("Partition the cis (intrachromosomal) BIOREP-level contacts into strata with these insert size boundaries "
               "and homogenize coverage within each stratum. Note that these are not the binning resolutions for the contact matrix. "
               "They are potential conjuncts used to homogenize the distribution of contacts during aggregation, prior to binning as a contact matrix.")
    )
    conditionReadConjuncts = param.ListSelector(objects = ["chrom1", "chrom2", "abs(pos1-pos2)", "cellID", "(rfrag1 == rfrag2)", "pair_type", "stratum"],
    check_on_set = False,
    doc = ("For statistics calling and coverage homogenization, Hich partitions the reads from each sample into blocks defined by the intersection of a set of "
           "CONJUNCTS, which are traits of the read. Reads having the same value for all CONJUNCTS will be assigned to the same READ BLOCK, and you can "
           "have Hich group your samples and downsample them to have the largest possible coverage while matching coverage across all READ BLOCKS. "
           "The benefit of this process is that it allows you to control for coverage differences across samples, as might result from different sequencing depths "
           "or techical errors during the assay. The conditionReadConjuncts parameter sets which CONJUNCTS should be used for this downsampling-based coverage control process "
           "specifically for condition-level samples."))
    conditionCisStrata = param.List(
        [100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000, 100000, 200000, 500000, 1000000, 2000000, 5000000],
        item_type = int,
        doc = ("Partition the cis (intrachromosomal) CONDITION-level contacts into strata with these insert size boundaries "
               "and homogenize coverage within each stratum. Note that these are not the binning resolutions for the contact matrix. "
               "They are potential conjuncts used to homogenize the distribution of contacts during aggregation, prior to binning as a contact matrix.")
    )
    techrepDownsampleToMeanDistribution = param.Boolean(False,
    doc = ("If set, computes the mean fraction of reads in each READ BLOCK defined by techrepReadConjuncts, ignoring outliers, and downsamples all techrep-level samples "
           "in each sample group (including outliers) to have a matching fraction of reads in each READ BLOCK."))
    biorepDownsampleToMeanDistribution = param.Boolean(False,
    doc = ("If set, computes the mean fraction of reads in each READ BLOCK defined by biorepReadConjuncts, ignoring outliers, and downsamples all biorep-level samples "
           "in each sample group (including outliers) to have a matching fraction of reads in each READ BLOCK."))
    conditionDownsampleToMeanDistribution = param.Boolean(False,
    doc = ("If set, computes the mean fraction of reads in each READ BLOCK defined by conditionReadConjuncts, ignoring outliers, and downsamples all condition-level samples "
           "in each sample group (including outliers) to have a matching fraction of reads in each READ BLOCK."))

    def to_dict(self):
        dict_repr = dict(self.param.values().items())
        dict_repr.pop("name")
        return dict_repr


aggregate_defaults = {
    "fullCoverageMerge": AggregateProfile(
        dedupSingleCell = False,
        dedupMaxMismatch = 3,
        dedupMethod = "max",
        techrepDedup = True,
        biorepDedup = True,
        conditionDedup = False,
        mergeTechrepToBiorep = True,
        mergeBiorepToCondition = True,
        techrepCisStrata = [],
        biorepCisStrata = [],
        conditionCisStrata = [],
        techrepReadConjuncts = [],
        biorepReadConjuncts = [],
        conditionReadConjuncts = [],
        techrepDownsampleToMeanDistribution = False,
        biorepDownsampleToMeanDistribution = False,
        conditionDownsampleToMeanDistribution = False
    ),
    "coverageMatchedMerge": AggregateProfile(
        dedupSingleCell = False,
        dedupMaxMismatch = 3,
        dedupMethod = "max",
        techrepDedup = True,
        biorepDedup = True,
        conditionDedup = False,
        mergeTechrepToBiorep = True,
        mergeBiorepToCondition = True,
        techrepDownsampleToMeanDistribution = True,
        biorepDownsampleToMeanDistribution = True,
        conditionDownsampleToMeanDistribution = True
    )
}

config["params"]["aggregate"].update(aggregate_defaults)

def view_aggregation_profiles():
    choices = to_nextflow_config(config["params"]["aggregate"]).split("\n") + ["Add new aggregation profile", "Main menu"]
    chosen = None
    action = None
    while chosen != "Main menu":
        chosen = Bullet(
                prompt = "These are the currently defined aggregation profiles. Select a profile name to edit or remove, or an individual parameter to view or modify.",
                choices = list(choices)
            ).launch()




def aggregation_menu():
    choices = {
        "View aggregation profiles.": view_aggregation_profiles,
        "Add aggregation profile.": None,
        "Remove aggregation profile.": None,
        "Main menu.": "main"
    }
    action = None

    while action != "main":
        chosen = Bullet(
            prompt = "Aggregation profiles are ways to downsample, control for sequencing depth and distance decay, deduplicate, and merge samples.",
            choices = list(choices.keys())
        ).launch()
        action = choices[chosen]
        if callable(action): action()

def tui_main_menu():
    choices = {
        "Default processing parameters for all samples.": None,
        "Special processing parameters for specific samples that override defaults.": None,
        "How to downsample, control for sequencing depth, deduplicate, and merge.": aggregation_menu,
        "Ways to select which samples are used for feature-calling.": None,
        "How to call compartment scores.": None,
        "How to call insulation scores.": None,
        "How to call loops.": None,
        "How to call A vs. B loop enrichments (diffloops).": None,
        "How to call HiC-Rep chromosome-wide correlations.": None,
        "When to generate read-level QC.": None,
        "How to run Hich appropriately in your computing environment.": None,
        "Where to save Hich outputs.": None,
        "Quit.": "quit"}
    action = None

    while action != "quit":
        chosen = Bullet(
            prompt = "Which aspect of your Hich workflow do you want to define?",
            choices = list(choices.keys())
        ).launch()
        action = choices[chosen]
        if callable(action): action()

"""
So for a config editing wizard, you can:

Define config subsections (which can be nested)
Config individual parameters
    Description associated with the parameter in general
    Description associated with the current parameter setting
    Optional?
    Use Param to define


Can we generalize over all subsections of config?
    We can display all the members of a config item
    Selecting an individual member lets you
        edit (but that's not as good as an editor),
        explain (but that's just referring to documentation),
    You could just walk them through the options one at a time
    Then give them the option to add it
    Also let them redo a profile (with current options selected by default)
    Or remove a profile
    Or rename a profile

The problem is it would be even harder to parse the textual config into classes that we can edit here.
We could potentially create a rigid structure in Python, then export it to Nextflow config but not have an import feature
"""

# from bullet import Bullet, Check, Input, YesNo, Numbers
import textwrap

from textual import on
from textual.app import App, ComposeResult
from textual.message import Message
from textual.widgets import Input as BaseInput, SelectionList, OptionList, Label
from textual.validation import ValidationResult, Validator
from textual.events import Key
import re

class RegexValidator(Validator):
    def __init__(self, pattern = None):
        self.pattern = pattern

    def validate(self, value: str) -> ValidationResult:
        """Check a string is equal to its reverse."""
        if self.pattern:
            success = bool(re.match(self.pattern, value))
            print(success)
        else:
            success = True
        return self.success() if success else self.failure(f"Does not conform to regex")

class Input(App):

    def __init__(self, prompt = None, default = "", indent = 0, strip = False, pattern = "", placeholder = "", *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.prompt = prompt
        self.default = default
        self.indent = indent
        self.strip = strip
        self.pattern = pattern
        self.placeholder = placeholder
        self.is_valid = None

    def compose(self) -> ComposeResult:
        if self.prompt is not None:
            self.label = Label(self.prompt)
            self.label.styles.width = 80
            yield self.label
    
        self.input = BaseInput(
            placeholder=self.placeholder or "",
            value = self.default,
            validators = [RegexValidator(self.pattern)]
        )
        yield self.input

    async def on_input_submitted(self, message: BaseInput.Submitted) -> None:
        # Capture the submitted value and exit the app
        submitted = message.value
        if self.strip:
            submitted = submitted.strip()
        if self.input.is_valid:
            self.exit(submitted)  # Exits the app and returns the value
    
    @on(BaseInput.Changed)
    def check_validation(self, event: BaseInput.Changed) -> None:
        self.is_valid = event.validation_result.is_valid

    def launch(self):
        return self.run()

class Check(App):
    def __init__(self, prompt = "", choices = [], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.prompt = prompt
        self.choices = choices
        self.label = None
        self.selection_list = None
    
    def compose(self) -> ComposeResult:
        if self.prompt is not None:
            self.label = Label(self.prompt)
            self.label.styles.width = 80
            yield self.label
    
        selections = [(choice, choice) for choice in self.choices]
        choice_type = type(self.choices[0])
        self.selection_list = SelectionList[choice_type](*selections)
        yield self.selection_list

    async def on_key(self, event: Key) -> None:
        # Check if the key pressed is "enter"
        if event.key == "enter":
            self.exit(self.selection_list.selected)

    def launch(self):
        return self.run()

class Bullet(App):
    def __init__(self, prompt = "", choices = [], *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.prompt = prompt
        self.choices = choices
        self.label = None
        self.option_list = None
    
    def compose(self) -> ComposeResult:
        if self.prompt is not None:
            self.label = Label(self.prompt)
            self.label.styles.width = 80
            yield self.label
    
        self.option_list = OptionList(*self.choices)
        yield self.option_list

    async def on_key(self, event: Key) -> None:
        # Check if the key pressed is "enter"
        if event.key == "enter":
            idx = self.option_list.highlighted
            choice = self.choices[idx]
            self.exit(choice)

    def launch(self):
        return self.run()

@dispatch(AggregateProfile)
def tui_create(profile):
    for param_name in profile.param:
        if param_name  == "name":
            continue

        param_value = profile.param.values()[param_name]
        param_obj = profile.param[param_name]
        doc = param_obj.doc #textwrap.fill(param_obj.doc, width = 70)
        prompt = f"{doc}\n{param_name}: "

        #print(doc)
        if param_obj.__class__ == param.Boolean:
            choices = ["Yes", "No"] if param_value else["No", "Yes"]
            choice = Bullet(prompt = prompt, choices = choices).launch()
            result = choice == "Yes"
        elif param_obj.__class__ == param.Integer:
            default = str(param_value)
            single_int = r"^\s*$|^\s*\d+\s*$"
            choice = None
            while choice is None:
                choice = Input(prompt = prompt, default = default, pattern = single_int).launch()
            result = param_value if not choice.strip() else int(choice)
        elif param_obj.__class__ == param.Selector:
            choices = param_obj.objects
            result = Bullet(prompt = prompt, choices = choices).launch()
        elif param_obj.__class__ == param.ListSelector:
            choices = list([str(x) for x in param_obj.objects])
            print(choices)
            result = Check(prompt = prompt, choices = choices).launch()
        elif param_obj.__class__ == param.List and param_obj.item_type == int:
            space_separated_int_list = r"^\s*$|^\d+( \d+)*$"
            default = " ".join([str(i) for i in param_value])
            # Create the Input object with the pattern
            int_list = None
            while int_list is None:
                int_list = Input(
                    prompt = prompt,
                    default = default,
                    pattern=space_separated_int_list
                ).launch()
            result = [int(i) for i in int_list.split()] if int_list.strip() else param_obj.default
        else:
            result = None
        if result is not None:
            setattr(profile, param_name, result)
    at_least_one_non_whitespace = r"\S+"
    name = Input(
        prompt = "Aggregate profile name: ",
        pattern = at_least_one_non_whitespace
    ).launch()
    
    print(name, end = "")
    print(to_nextflow_config(profile))


tui_create(AggregateProfile())