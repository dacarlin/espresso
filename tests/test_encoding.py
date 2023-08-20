import coaster 


protein_1 = "MENFHHRPFKGGFGVGRVPTSLYYSLSDFSLSAISIFPTHYDQPYLNEAPSWYKYSLES"
protein_2 = "MACDEFGHIKLMNPQRSTVWYMACDEFGHIKLMNPQRSTVWY"


def test_encode_sequence_of_met():
    protein = "MMM"
    expected = "ATGATGATG"
    assert coaster.encode_sequence(protein, "sc") == expected 


def test_different_built_in_tables():
    tables = ["ec", "sc", "yl"]
    for table in tables:
        seq = coaster.encode_sequence(protein_1, table)
        assert seq[:3] == "ATG"
        seq = coaster.encode_sequence(protein_2, table)
        assert len(seq) / 3 == len(protein_2)

