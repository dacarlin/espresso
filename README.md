# Espresso

Use generative AI to create highly-expressed synthetic genes


## User guide 

### Installation 

For now, install Espresso by cloning the repository and using your favorite
Python package management system. For example, using a regular virtual
environment, first create and activate your environment (here I am using 
the Fish shell)  

```fish 
python -m venv .venv 


```shell 
git clone https://github.com/dacarlin/espresso.git
cd espresso
python -m pip install -e .  
```


### Quickstart 

First, make sure you can import Espresso! 

```python 
import espresso 
```

Next, let's try generating a coding sequence for _E. coli_ using a protein sequence 
that we provide 

```python
my_protein = "MENFHHRPFKGGFGVGRVPTSLYYSLSDFSLSAISIFPTHYDQPYLNEAPSWYKYSLESGLVCLYLYLIYRWITRSF"
my_gene = espresso.design_coding_sequence(my_protein, model="ec")
```

The resulting sequence will be optimized for expression in _E. coli_. Read on to learn 
more about Espresso and how to perform more complicated tasks. 


### Data for the independent model

For the independent model, learn the codon frequency from a set of genes 
using the provided `learn.py` script. For example 

```bash 
python learn.py espresso/data/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz	
```

This will output a JSON file that can be used with the `IndependentEncoder` class


### Data for the transformer models 

For the transformer models, learn whole sequence design from native sequences 
stored in a FASTA file 
