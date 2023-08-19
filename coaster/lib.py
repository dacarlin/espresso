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
