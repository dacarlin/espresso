import espresso 


protein = "MENFHHRPFKGGFGVGRVPTSLYYSLSDFSLSAISIFPTHYDQPYLNEAPSWYKYSLESGLV"


def test_encoding_interface():
    enc = espresso.design_coding_sequence(protein, model="sc")
    assert enc.startswith("ATG")


def test_scrubbing_interface():

    # N, codons AAU AAC only 
    # the peptide KN can be encoded 4 ways 
    # 1. AAA AAT 
    # 2. AAA AAC
    # 3. AAG AAT
    # 4. AAG AAC 
    # so for example if you ban a 5-A, then you only get 3 and 4, 
    # and if you ban 5-A and GAAC, you will always get AAG AAT (3) 

    avoid = [
        "AAAAA", 
        "GAAC", 
    ]
    enc = espresso.scrub_sequence("AAAAAT", avoid=avoid, model="sc")
    assert enc == "AAGAAT"

