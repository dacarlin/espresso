# coaster

Use generative AI to create highly-expressed synthetic genes


## User guide 

### Installation 

### Quickstart 

```python 
import coaster 


protein = "MENFHHRPFKGGFGVGRVPTSLYYSLSDFSLSAISIFPTHYDQPYLNEAPSWYKYSLESGLVCLYLYLIYRWITRSF"
gene = coaster.encode_sequence(protein) 
```

### Data for the independent model

For the independent model, learn the codon frequency from a set of genes 
using the provided `learn.py` script. For example 

```bash 
python learn.py coaster/data/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz	
```

This will output a JSON file that can be used with the `IndependentEncoder` class


### Data for the transformer models 

For the transformer models, learn whole sequence design from native sequences 
stored in a FASTA file 
