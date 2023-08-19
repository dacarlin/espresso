from coaster.lib import IndependentEncoder 
from coaster.data import sc_codon_use, ec_codon_use, yl_codon_use


def encode_sequence(protein_sequence, model="sc"):
    """Create a gene sequence from a protein sequence

    Parameters
    ----------
    protein_sequence: str
        Protein sequence as a string 
    codon_table: str
        The name of a codon model in Coaster

    Examples
    --------
    An example of encoding a protein 

    >>> encoded = encode_sequence("MMM")
    """

    choices = {
        "sc": IndependentEncoder(sc_codon_use), 
        "ec": IndependentEncoder(ec_codon_use),
        "yl": IndependentEncoder(yl_codon_use),
    }

    if model in choices:
        model = choices[model] 
    else:
        raise ValueError(f'Model "{model}" not found')

    sequence = model.generate_sequence(protein_sequence) 

    return sequence 
    

def scrub_sequence(nucleotide_sequence, codon_table="sc"):
    return nucleotide_sequence 


