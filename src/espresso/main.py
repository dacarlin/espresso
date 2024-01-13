from espresso.lib import IndependentEncoder, Scrubber, AvoidMotif
from espresso.data import sc_codon_use, ec_codon_use, yl_codon_use


def design_coding_sequence(protein_sequence, model="sc"):
    """Create a gene sequence from a protein sequence

    Parameters
    ----------
    protein_sequence: str
        Protein sequence as a string 
    codon_table: str
        The name of a codon model

    Examples
    --------
    An example of encoding a protein 

    >>> encoded = design_coding_sequence("MMM")
    """

    choices = {
        "sc": IndependentEncoder(sc_codon_use), 
        "ec": IndependentEncoder(ec_codon_use),
        "yl": IndependentEncoder(yl_codon_use),
        #"fungi-v1": TransformerEncoder(fungi_v1), 
    }

    if model in choices:
        model = choices[model] 
    else:
        raise ValueError(f'Model "{model}" not found')

    sequence = model.generate_sequence(protein_sequence) 

    return sequence 
    

def scrub_sequence(nucleotide_sequence, avoid, model="sc"):
    """Scrub a nucleotide sequence of specific motifs or regions of GC content"""

    # process avoid options 
    avoid = [AvoidMotif(x) for x in avoid]

    # process model options 
    choices = {
        "sc": IndependentEncoder(sc_codon_use), 
        "ec": IndependentEncoder(ec_codon_use),
        "yl": IndependentEncoder(yl_codon_use),
    }

    if model in choices:
        model = choices[model] 
    else:
        raise ValueError(f'Model "{model}" not found')

    # scrub 
    scrubber = Scrubber(avoid=avoid, model=model)

    return scrubber.scrub(nucleotide_sequence) 

