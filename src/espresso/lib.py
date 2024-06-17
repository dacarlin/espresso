import re 
from itertools import product

import torch 
import numpy
from biotite.sequence import NucleotideSequence 
import torchtext; torchtext.disable_torchtext_deprecation_warning()
from torchtext.vocab import build_vocab_from_iterator

from espresso.model import Seq2SeqTransformer


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


class IndependentModel:
    """Uses codon usage data to create encodings 
    with the same global codon frequency as the 
    input data"""

    def __init__(self, codon_use_data, threshold=0.):
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

        # remove entries below threshold 
        self.table[self.table < threshold] = 0.

        # re-normalize the mapping table by the row 
        self.table = self.table / self.table.sum(axis=1, keepdims=True)

    def generate_sequence(self, protein_sequence):
        """Generate a nucletide coding sequence for a provided protein sequence"""

        sequence = ""
        for residue in protein_sequence:
            idx = RESIDUE_TO_INDEX[residue]
            p = self.table[idx]
            codon = numpy.random.choice(CODONS, p=p)
            sequence += codon 

        return sequence 
    

class TopCodonModel:
    """Uses codon usage data to create encodings 
    with the single most common codon for each position"""

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

        sequence = ""
        for residue in protein_sequence:
            idx = RESIDUE_TO_INDEX[residue]
            p = self.table[idx]
            codon = CODONS[numpy.argmax(p)]
            sequence += codon 

        return sequence 


class TransformerModel:
    """Uses a pre-trained transformer model to design coding sequences"""
    def __init__(self, model_path):

        protein_vocab = list("ACDEFGHIKLMNPQRSTVWY*")
        codon_vocab = list("".join(x) for x in product("ATCG", repeat=3))

        SRC_LANGUAGE = "protein"
        TGT_LANGUAGE = "codons"

        # Place-holders
        token_transform = {}
        vocab_transform = {}

        UNK_IDX, PAD_IDX, BOS_IDX, EOS_IDX = 0, 1, 2, 3
        special_symbols = ['<unk>', '<pad>', '<bos>', '<eos>']
        self.BOS_IDX = BOS_IDX
        self.EOS_IDX = EOS_IDX

        # first the protein vocab!
        vocab_transform[SRC_LANGUAGE] = build_vocab_from_iterator([protein_vocab], min_freq=1, specials=special_symbols, special_first=True)

        # now the codon vocab!
        vocab_transform[TGT_LANGUAGE] = build_vocab_from_iterator([codon_vocab], min_freq=1, specials=special_symbols, special_first=True)

        for lang in [SRC_LANGUAGE, TGT_LANGUAGE]:
            vocab_transform[lang].set_default_index(UNK_IDX)

        # now the token transforms
        token_transform[SRC_LANGUAGE] = lambda x: x.split(" ")
        token_transform[TGT_LANGUAGE] = lambda x: x.split(" ")

        self.vocab_transform = vocab_transform
        self.token_transform = token_transform

        # Create a new empty instance of the model 
        # 
        # In our training, detailed in `espresso/learn`, we chose 
        # the following model params:
        #
        # encoder layers == decoder layers, 3
        # model dim, 64 
        self.model = Seq2SeqTransformer(3, 3, 64, 8, len(vocab_transform[SRC_LANGUAGE]), len(vocab_transform[TGT_LANGUAGE]), 256) 

        # load the specified trained model 
        state_dict = torch.load(model_path, map_location=torch.device('cpu'))
        self.model.load_state_dict(state_dict)

        # set in eval mode 
        self.model.eval() 

    def sequential_transforms(self, *transforms):
        def func(txt_input):
            for transform in transforms:
                txt_input = transform(txt_input)
            return txt_input
        return func

    # function to add BOS/EOS and create tensor for input sequence indices
    def tensor_transform(self, token_ids):
        return torch.cat((torch.tensor([self.BOS_IDX]),
                        torch.tensor(token_ids),
                        torch.tensor([self.EOS_IDX])))

    def generate_sequence(self, protein_sequence, verbose=False):
        """Generate a CDS for provided protein sequence using a generative model
        
        Raises
        ------
        AssertionError
            If, after 5 samples, none of the generated sequences translate to 
            the provided protein sequence
        """

        text_transform = {}
        for ln in ["protein", "codons"]:
            text_transform[ln] = self.sequential_transforms(self.token_transform[ln], 
                                                    self.vocab_transform[ln], 
                                                    self.tensor_transform) 
        src_sentence = " ".join(protein_sequence) 
        src = text_transform["protein"](src_sentence).view(-1, 1)
        num_tokens = src.shape[0]
        src_mask = (torch.zeros(num_tokens, num_tokens)).type(torch.bool)

        n_iter = 0
        max_iter = 20
        while n_iter < max_iter:
            n_iter += 1
            tgt_tokens = self.model.sample(src, src_mask, max_len=num_tokens + 1, start_symbol=self.BOS_IDX).flatten() 
            result = " ".join(self.vocab_transform["codons"].lookup_tokens(list(tgt_tokens.cpu().numpy()))).replace("<bos>", "").replace("<eos>", "")
            sequence = NucleotideSequence(result.replace(" ", ""))
            trimmed_sequence = sequence[:num_tokens * 3]
            translation = trimmed_sequence.translate(complete=True)
            if str(translation) == str(protein_sequence):
                break  # break out of the `while` loop 

        if verbose:
            print("espresso: input sequence length including <bos> and <eos>:", len(src))
            print("espresso: raw codons looked up from tokens", result)
            print("espresso: joined sequence", sequence) 
            print("espresso: joined sequence trimmed to length", trimmed_sequence) 
            print("espresso: translation of trimmed sequence", translation) 

        if not str(translation) == str(protein_sequence):
            raise RuntimeError(f"espresso: The model produced a nucleotide sequence that doesn't translate back to the original nucleotide sequence after {max_iter} iterations")

        return str(trimmed_sequence) 
        

class Scrubber:
    """Uses a codon model to scrub a sequence of undesired characteristics
    like specific motifs, GC content, etc"""

    def __init__(self, avoid, model):
        self.avoid = avoid 
        self.model = model 

    def scrub(self, nucleotide_sequence, max_iterations=50):
        """For `max_iterations`, attempt to generate a new sequence
        without any of the undesired features"""
       
        iterations = 0 
        for n in range(max_iterations):
            nucleotide_sequence = self.resample_sequence(nucleotide_sequence)
            if not self.identify_codons_to_resample(nucleotide_sequence):
                return nucleotide_sequence 
            iterations += 1 
            if iterations == 50:
                raise RuntimeError((f"Reached maximum allowed iterations {max_iterations} when"
                                     "encoding a sequence starting with {nucleotide_sequence[:24]}"
                                     "with {len(self.avoid)} constraints"))

    def identify_codons_to_resample(self, nucleotide_sequence):
        bases = []
        for thing in self.avoid:
            my_bases = thing(nucleotide_sequence) 
            bases.extend(my_bases) 
        bases = list(set(bases)) 
        codons_to_resample = list(set(x // 3 for x in bases))
        return codons_to_resample

    def resample_sequence(self, nucleotide_sequence):
        # create a sequence object to get the original codons 
        sequence = CodonSequence(nucleotide_sequence) 

        # score the sequence with each avoid block 
        codons_to_resample = self.identify_codons_to_resample(nucleotide_sequence) 
       
        # resample the specified codons (or generate a new sequence from the model!) 
        # probably could rethink this as a "mask" for the model 
        new_sequence = ""
        for idx in range(len(sequence.codons)):
            if idx in codons_to_resample:
                residue = translate(sequence.codons[idx])
                new_codon = self.model.generate_sequence(residue)
                new_sequence += new_codon
            else:
                new_sequence += sequence.codons[idx] 

        return new_sequence 


class AvoidMotif:
    """Avoid a particular sequence
    
    Examples
    --------
    To get the positions (base pairs) that are involved in matching
    the motif you wanna avoid 
    >>> sequence = "ATGCCCCCC"
    >>> motif = AvoidMotif("CCCCCC")
    >>> motif(sequence)  
    # [3, 4, 5, 6, 7, 8] 
    """
    def __init__(self, motif):
        self.motif = motif 
    
    def __call__(self, sequence):
        # get a list of the positions involved 
        positions = []
        for m in re.finditer(self.motif, sequence):
            my_range = list(range(m.start(), m.start() + len(self.motif)))
            positions.extend(my_range)

        return list(set(positions))
