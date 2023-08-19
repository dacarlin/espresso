import coaster 


def test_encode_sequence_of_met():
    protein = "MMM"
    expected = "ATGATGATG"
    assert coaster.encode_sequence(protein, "sc") == expected 
