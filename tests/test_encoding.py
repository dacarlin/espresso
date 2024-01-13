import espresso 
from espresso.data import ec_codon_use, fungi_v1
from espresso.lib import TopCodonEncoder, IndependentEncoder, TransformerEncoder


protein_1 = "MENFHHRPFKGGFGVGRVPTSLYYSLSDFSLSAISIFPTHYDQPYLNEAPSWYKYSLES"
protein_2 = "MACDEFGHIKLMNPQRSTVWYMACDEFGHIKLMNPQRSTVWY"


def test_top_codon_encoder():
    model = TopCodonEncoder(ec_codon_use)
    seq = model.generate_sequence(protein_1)
    assert seq[:3] == "ATG"

def test_top_codon_encoder():
    model = IndependentEncoder(ec_codon_use)
    seq = model.generate_sequence(protein_1)
    assert seq[:3] == "ATG"


def test_design_coding_sequence_of_met():
    protein = "MMM"
    expected = "ATGATGATG"
    assert espresso.design_coding_sequence(protein, "sc") == expected 


def test_different_built_in_independent_tables():
    tables = ["ec", "sc", "yl"]
    for table in tables:
        seq = espresso.design_coding_sequence(protein_1, table)
        assert seq[:3] == "ATG"
        seq = espresso.design_coding_sequence(protein_2, table)
        assert len(seq) / 3 == len(protein_2)


def test_transformer_encoder():
    model = TransformerEncoder(fungi_v1)
    seq = model.generate_sequence(protein_2)
    assert seq[:3] == "ATG"