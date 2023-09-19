import espresso 


protein_1 = "MENFHHRPFKGGFGVGRVPTSLYYSLSDFSLSAISIFPTHYDQPYLNEAPSWYKYSLES"
protein_2 = "MACDEFGHIKLMNPQRSTVWYMACDEFGHIKLMNPQRSTVWY"


def test_encode_sequence_of_met():
    protein = "MMM"
    expected = "ATGATGATG"
    assert espresso.encode_sequence(protein, "sc") == expected 


def test_different_built_in_tables():
    tables = ["ec", "sc", "yl"]
    for table in tables:
        seq = espresso.encode_sequence(protein_1, table)
        assert seq[:3] == "ATG"
        seq = espresso.encode_sequence(protein_2, table)
        assert len(seq) / 3 == len(protein_2)

