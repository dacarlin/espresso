import io 

from espresso.lib import TopCodonModel, IndependentModel, TransformerModel, Scrubber, AvoidMotif
from espresso.data import sc_codon_use, ec_codon_use, yl_codon_use, fungi_v1


def get_choices():
    """Get a dict of all available models """
    # choices = {
    #     "yeast-top": TopCodonModel(sc_codon_use), 
    #     "coli-top": TopCodonModel(ec_codon_use), 
    #     "sc": IndependentModel(sc_codon_use), 
    #     "ec": IndependentModel(ec_codon_use),
    #     "yl": IndependentModel(yl_codon_use),
    #     "fungi-v1": TransformerModel(fungi_v1), 
    # }
    choices = {
        "yeast-top": (TopCodonModel, sc_codon_use), 
        "coli-top": (TopCodonModel, ec_codon_use), 
        "sc": (IndependentModel, sc_codon_use), 
        "ec": (IndependentModel, ec_codon_use),
        "yl": (IndependentModel, yl_codon_use),
        "fungi-v1": (TransformerModel, fungi_v1), 
    }
    return choices


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

    choices = get_choices() 

    if model in choices.keys():
        model_cls, model_data = choices[model] 
        model = model_cls(model_data)
        
        # For the transformer models, the data is stored on disk as 
        # a Torch tensor, and read into a `BytesIO` wrapper to be used.
        # Because of this, we need to seek the `BytesIO` object to the
        # beginning each time 
        if isinstance(model_data, io.BytesIO):
            model_data.seek(0)
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
        "sc": IndependentModel(sc_codon_use), 
        "ec": IndependentModel(ec_codon_use),
        "yl": IndependentModel(yl_codon_use),
    }

    if model in choices:
        model = choices[model] 
    else:
        raise ValueError(f'Model "{model}" not found')

    # scrub 
    scrubber = Scrubber(avoid=avoid, model=model)

    return scrubber.scrub(nucleotide_sequence) 

