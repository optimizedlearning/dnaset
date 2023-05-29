# Utilities for Loading Genomic Data in Pytorch

on the scc, you will need to `module load python3/3.10.5 pytorch ucscutils bedtools gcc`
and then set up a virtual environment (or `source`) a previous one:
```
python -m venv dnaset_env
source dnaset_env/bin/activate
pip install --upgrade  pip
pip install -r requirements.txt
```

Note that we need a more recent python version, otherwise some of the huggingface libraries will complain with some SSL errors.

You can do all of this with a `source scc_setup.sh`. Note that the first time this runs it will take a long time (because the `pip install`
will compile a bunch of comp bio tools).

Trained tokenizer files are located in `/projectnb/aclab/dnaset/saved_tokenizers/`. For example:
```
from transformers import PreTrainedTokenizerFast
tokenizer = PreTrainedTokenizerFast.from_pretraine('/projectnb/aclab/dnaset/saved_tokenizers/multispecies.32768')
tokenizer.encode('ACCTCAGCATTGACTCTA') # -> [11530, 17266, 161]
```