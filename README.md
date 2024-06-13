# Espresso

Espresso is a free, open-source app that allows anyone to create coding sequences 
for synthetic genes. Espresso supports a variety of different host organisms and 
mRNA design algorithms. Using Espresso, you can harness the power of generative AI 
to create highly-expressed synthetic genes for your host of choice. 

- **New in Espresso 1.1: pre-trained transformer models for yeast!**
    - Use pre-trained transformer models to design genes for expression in yeast 
    - Simply specify the model `fungi-v1` when designing genes 
    - Example usage: `import espresso; espresso.design_coding_sequence("MSE…NST", "fungi-v1")` to use a pre-trained transformer model to design a synthetic gene for the provided sequence 


## User guide 

### Installation 

Install Espresso by cloning the repository and using your favorite
Python package management system. For example, using a regular virtual
environment, first create and activate your environment (here I am using 
the Fish shell)  

```shell 
# create and activate virtual environment 
python -m venv .venv 
source .venv/bin/activate.fish
```

and then clone and install with Pip

```shell 
# clone and install with pip 
git clone https://github.com/dacarlin/espresso.git
cd espresso
python -m pip install -e .  
```

**Coming soon (summer 2024):** Espresso will soon be available as a web app. 


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

### Designing coding sequences 

There are many important parts of designing a sequence for expression in a specified host. We must specify the sequence, the host organism, and model we want to use, and also any constraints on the problem—for example, avoiding specific restriction enzyme sites. 

Espresso provides easy access to common workflows, while maintaining the flexibility to specify complex constraints. 

#### Specifying input sequences 

Specify your protein sequence as a string. 

```python 
my_protein = "MENFHHRPFKGGFGVGRVPTSLYYSLSDFSLSAISIFPTHYDQPYLNEAPSWYKYSLESGLVCLYLYLIYRWITRSF"
```

#### Choosing a design model 

Available models include two types: there are "independent" models, which produce genes which recapitulate the observed codon frequencies and where sites are independent. And there are "transformer" models, which design genes using a pre-trained language model that models the distribution of natural sequences. 

| Model slug     | Model type | Description      |
| --- | --- | ---- | 
| sc    | Codon frequency with threshold | Independent model trained on _Saccharomyces cerevisiae_ native genes  |
| ec    | Codon frequency with threshold | Independent model trained on _Escherichia coli_ (_E. coli_) native genes  |
| sc    | Codon frequency with threshold | Independent model trained on _Yarrowia lipolytica_ native genes  |
| yeast-top | Top codon only | Model trained on _Saccharomyces cerevisiae_ native genes that outputs deterministic sequences containing only the most frequent codon for each residue |
| coli-top | Top codon only | Model trained on _Eschericia coli_ (_E. coli_) native genes that outputs deterministic sequences containing only the most frequent codon for each residue |
| fungi-v1 | Transformer | Deep neural network with a transformer architecture trained on native genes from thousands of fungal taxa |

To use the design model when designing a sequence, pass the slug of the model to `design_coding_sequence`. For example: 

```python 
# design for E. coli 
espresso.design_coding_sequence(my_protein, "ec")

# design using transformer model for yeast 
espresso.design_coding_sequence(my_protein, "fungi-v1")
```

#### Scrubbing sequences of undesired motifs 

Often, we want to avoid specific motifs in our designed sequences. To solve this problem, Espresso implements a "scrubber", which "scrubs" sequences of specific motifs in a manner similar to image inpainting. 

The algorithm is simple: Espresso scans a provided sequence for specified motifs that should be avoided, and then inpaints a new sequence using a host-specific model. I designed this algorithm to be a fast and accurate solution to avoiding restriction enzyme sites while maintaining optimal biological coding. It's also quite fast. 

```python 
avoid = [
    "AAAAAA",  # avoid poly-A, synthesis restriction  
    "GAATTC",  # avoid EcoRI and its reverse complement  
]

protein = "MENFHHRPFKGGFGVGRVPTSLYYSLSDFSLSAISIFPTHYDQP"

# first design a candidate 
candidate = espresso.design_coding_sequence(protein, model="ec")

# then scrub the candidate 
scrubbed = espresso.scrub_sequence(candidate, avoid, model="ec")
```


### Training your own models 

#### Training your own codon models 

For the independent model, learn the codon frequency from a set of genes 
using the provided `count_codons.py` script. For example 

```bash 
python learn.py espresso/data/cds/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa.gz	
```

This will output a JSON file that can be used with the `IndependentEncoder` class. Save the JSON file to disk, and then create your own encoder class.

```python
from espresso.lib import IndependentModel

path_to_json = "path/to/my_codons.json"
with open(path_to_json) as fn:
    data = json.read(fn)

my_model = IndependentModel(data)
```

You can now use the model's functions with your custom dataset, such as 

```python
my_model.generate_sequence("MSENT")
```


#### Details on training the transformer models 

For more details on how the transformer models are trained, along with the code implementation, please see [my blog post](https://alexcarlin.bearblog.dev/using-generative-ml-to-design-native-looking-genes-new/).
