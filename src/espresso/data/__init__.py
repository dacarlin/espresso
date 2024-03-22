import io 
import json 
import pkgutil


sc_codon_use = json.loads(pkgutil.get_data("espresso", "data/independent/sc.json"))
ec_codon_use = json.loads(pkgutil.get_data("espresso", "data/independent/ec.json"))
yl_codon_use = json.loads(pkgutil.get_data("espresso", "data/independent/yl.json"))
fungi_v1 = io.BytesIO(pkgutil.get_data("espresso", "data/transformer/fungi-v1.pt"))