from biotite.sequence import NucleotideSequence 
from itertools import product

import numpy 


STOP_CODONS = ["TAA", "TAG", "TAA", "TGA"]
RESIDUE_TO_INDEX = dict(zip("ACDEFGHIKLMNPQRSTVWY", range(20)))
INDEX_TO_RESIDUE = {v: k for k, v in RESIDUE_TO_INDEX.items()}
CODON_TO_INDEX = dict(zip(list("".join(x) for x in product("ATCG", repeat=3)), range(64)))
INDEX_TO_CODON = {v: k for k, v in CODON_TO_INDEX.items()}
CODONS = list(CODON_TO_INDEX.keys())


def translate(nucleotide_sequence):
    """Helper function to translate codons into residues"""
    dna = NucleotideSequence(nucleotide_sequence)
    return str(dna.translate(complete=True))


class CodonSequence:
    """A nucleotide sequence consisting of one or more codons"""

    def __init__(self, sequence):
        """Create a new CodonSequence after verifying data"""
        if not len(sequence) % 3 == 0:
            raise ValueError(f"Sequence length must be a multiple of 3, not {len(sequence)}")

        self.codons = list(sequence[x:x + 3] for x in range(0, len(sequence), 3))
        if len(self.codons) < 1:
            raise ValueError("A sequence must have at least 3 bp") 

        self.sequence = "".join(self.codons) 


class CodonTable:
    """A codon table for a particular task, or organism"""

    def __init__(self, frame):
        assert "codon" in frame.columns, "Frame must have a 'codon' column" 
        assert "count" in frame.columns, "Frame must have a 'count' column"

        self.codons = list(frame["codon"])
        self.counts = list(frame["count"])


class IndependentEncoder:
    """Uses codon usage data to create encodings 
    with the same global codon frequency as the 
    input data"""

    def __init__(self, codon_use_data):
        # remove stop codons
        self.codon_use_data = {k: v for k, v in codon_use_data.items() if k not in STOP_CODONS}

        # construct mapping table 
        self.table = numpy.zeros((20, 64))
        for codon, count in self.codon_use_data.items():
            my_codon_index = CODON_TO_INDEX[codon]
            my_residue_index = RESIDUE_TO_INDEX[translate(codon)]
            self.table[my_residue_index, my_codon_index] = count 

        # normalize the mapping table by the row 
        self.table = self.table / self.table.sum(axis=1, keepdims=True)

    def generate_sequence(self, protein_sequence):
        """Generate a nucletide coding sequence for a provided protein sequence"""

        # later, you can use the mapping table to get a probability dist over the 64 possible codons 
        # get your residue index, idx 
        # index into the table, p = table[idx]
        # p is a 64 probability distribution over the amino acids 
        # pick one, get an index, 32 
        # get the codon INDEX_TO_CODON[32]

        sequence = ""
        for residue in protein_sequence:
            idx = RESIDUE_TO_INDEX[residue]
            p = self.table[idx]
            codon = numpy.random.choice(CODONS, p=p)
            sequence += codon 

        return sequence 


        
