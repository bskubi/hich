from collections import Counter
from hich.parse.pairs_file import PairsFile
from numpy import searchsorted
from dataclasses import dataclass, field
from hich.stats.selection_sampler import SelectionSampler
from hich.stats.pairs_rv import PairsRV
import sys

@dataclass
class PairsSampler(SelectionSampler):
    rv: PairsRV = None

    def parse_events(self, pairs_file: PairsFile, max_pairs = float('inf')):
        for i, pair in enumerate(pairs_file):
            if i >= max_pairs: break

            event = rv.event(pair)
            yield pair, event

    def count_events(self, pairs_file: PairsFile, max_pairs = float('inf')):
        for pair, event in self.parse_events(pairs_file, max_pairs):
            self.count(event)
            yield event
    
    def sample_events(self, pairs_file: PairsFile, max_pairs = float('inf')):
        for pair, event in self.parse_events(pairs_file, max_pairs):
            keep = self.sample(event)
            if keep: yield pair, event
    
    def write_sample(self, input_pairs: PairsFile, output_filename: str, max_pairs = float('inf')):
        header = input_pairs.header
        command_to_generate_pairs_file = ' '.join(sys.argv)
        header.command.append(command_to_generate_pairs_file)
        output_pairs = PairsFile(output_filename, mode = "w", header = header)
        
        for pair, event in self.sample_events(input_pairs, max_pairs):
            output_pairs.write(pair)
