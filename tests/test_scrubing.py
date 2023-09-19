from espresso.lib import Scrubber, AvoidMotif, IndependentEncoder
from espresso.data import sc_codon_use 


def test_avoid_motif():
    motif = AvoidMotif("AAAAA")
    positions = motif("AAAAAGAAAAGG")
    assert positions == [0, 1, 2, 3, 4] 

    positions = motif("GGGAAAAACTCTCTCATAT")
    assert positions == [3, 4, 5, 6, 7] 


def avoid_two_motifs():
    motif_1 = AvoidMotif("AAAAA")
    motif_2 = AvoidMotif("GAAC")
    
    motifs = [
        motif_1, 
        motif_2, 
    ]

    sequence = "AAAAAT"
    assert motif_1(sequence) == [0, 1, 2, 3, 4]

    assert motif_2(sequence) == []


def test_scrubber_object():
    # K, codons AAA AAG only 
    # N, codons AAU AAC only 
    # the peptide KN can be encoded 4 ways 
    # 1. AAA AAT 
    # 2. AAA AAC
    # 3. AAG AAT
    # 4. AAG AAC 
    # so for example if you ban a 5-A, then you only get 3 and 4, 
    # and if you ban 5-A and GAAC, you will always get AAG AAT (3) 

    avoid = [
        AvoidMotif("AAAAA"), 
        AvoidMotif("GAAC"), 
    ]
    model = IndependentEncoder(sc_codon_use) 
    scrubber = Scrubber(avoid=avoid, model=model)
    result = scrubber.scrub("AAAAAT")  # dipeptide KN
    assert result == "AAGAAT"
