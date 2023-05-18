

from typing import List, Tuple, Optional, Any
import itertools
import util

import torch
import transformers
from transformers import PreTrainedTokenizer
from tokenizers import (
        decoders,
        models,
        normalizers,
        pre_tokenizers,
        processors,
        trainers,
        Tokenizer
        )
from pybedtools import BedTool
from pyfaidx import Fasta


def batch_iterator(it, batch_size):
    '''
    converts an iterator into one that returns "batches"
    of size batch_size.    
    '''
    while(True):
        batch = list(itertools.islice(it, 0, batch_size))
        batch = [str(s.sequence) for s in batch]
        yield batch
        if len(batch)< batch_size:
            break



def tokenizer_from_bed(
        bed_files: List[str],
        fasta_files = List[str],
        vocab_size=25000,
        ):
    '''
    generates a tokenizer from bed files.

    Args:
        bed_files: a list of bed files specifying the sequences
            to use to train the tokenizer.
        fasta_files: a list of fasta files corresponding to
            bed_files. fasta_files[i] is the reference for
            bed_files[i].

    returns:
        a Tokenizer object containing a BPE tokenizer.
    '''

    sequence_generators = [
            util.bed_to_sequence_generator(bed, fasta) for  bed,fasta in zip(bed_files, fasta_files)
            ]

    generator = batch_iterator(
            itertools.chain(*sequence_generators),
            batch_size=1000
            )


    tokenizer = Tokenizer(models.BPE(unk_token="[UNK]"))

    tokenizer.normalizer = normalizers.Lowercase()

    special_tokens = ["[UNK]", "[PAD]", "[CLS]", "[SEP]", "[MASK]"]
    trainer = trainers.BpeTrainer(vocab_size=vocab_size, special_tokens=special_tokens, show_progress=True)

    tokenizer.train_from_iterator(generator, trainer=trainer)

    return tokenizer



