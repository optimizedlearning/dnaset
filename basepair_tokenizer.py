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
from datasets import Dataset

import random

def batch_iterator(dataset, batch_size, num_bed_lines):
    '''
    converts an iterator into one that returns "batches"
    of size batch_size.    
    '''
    for i in range(0, len(dataset), int(len(dataset) *  batch_size/num_bed_lines)):
        
        yield dataset[i:i+batch_size]['sequence']


def build_tokenizer(
        fasta_files: List[str],
        bed_files: List[str] = [],
        vocab_size=25000,
        fasta_chunk_size=100000,
        num_bed_lines=1000000,
        subsample=None
        ):
    '''
    generates a tokenizer from bed files.

    Args:
        fasta_files: a list of fasta files to tokenizer
        bed_files: a list of bed files specifying the sequences
            to use to train the tokenizer. bed_files[i]
            uses fasta_files[i] as reference.
            If len(fasta_files)>len(bed_files), then any
            remaining fasta files are tokenizer completely.

    returns:
        a Tokenizer object containing a BPE tokenizer.
    '''

    def generator():
        sequence_generators = [
                util.bed_to_sequence_generator(bed, fasta, True) for  bed,fasta in zip(bed_files, fasta_files[:len(bed_files)])
                ]
        sequence_generators.extend(
            util.fasta_chunk_generator(fasta, fasta_chunk_size, True) for fasta in fasta_files[len(bed_files):]
        )

        return itertools.chain(*sequence_generators)

    dataset = Dataset.from_generator(generator, cache_dir='temp_cache/')
    if subsample is not None:
        num_bed_lines = int(len(dataset) * subsample)
    generator = batch_iterator(
            dataset,
            batch_size=100,
            num_bed_lines=num_bed_lines,
            )


    tokenizer = Tokenizer(models.BPE(unk_token="[UNK]"))

    tokenizer.normalizer = normalizers.Lowercase()
    tokenizer.pre_tokenizer = pre_tokenizers.ByteLevel(False,False)
    
    tokenizer.decoder = decoders.ByteLevel()
    
    special_tokens = ["[UNK]", "[PAD]", "[CLS]", "[SEP]", "[MASK]"]
    trainer = trainers.BpeTrainer(vocab_size=vocab_size, special_tokens=special_tokens, show_progress=True)

    tokenizer.train_from_iterator(generator, trainer=trainer)

    
    return tokenizer

    


