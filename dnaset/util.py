import numpy as np
import torch

from pybedtools import BedTool
from pyfaidx import Fasta

def seq_to_array(s, copy=True):
    array = np.asarray(memoryview(str(s).encode('ascii')), dtype=np.uint8)
    if copy:
        array = array.copy()
    return array

def seq_to_torch(s, copy=True):
    array = torch.Tensor(memoryview(str(s).encode('ascii')), dtype=torch.uint8)
    if copy:
        array = array.copy()
    return array


def char_view(s):
    return s.view(dtype='|S1')

def int_view(s):
    return s.view(dtype=np.uint8)

def fasta_chunk_generator(
    fasta: str,
    chunk_size: int,
    seq_as_str: bool = False,
):
    fasta = Fasta(fasta)
    chroms = fasta.keys()
    for chrom in chroms:
        for start in range(0, len(fasta[chrom]), chunk_size):
            stop = min(len(fasta[chrom]), start+chunk_size)
            sequence = fasta[chrom][start:stop]
            if seq_as_str:
                sequence = str(sequence)
            yield {
                'chrom': chrom,
                'start': start,
                'stop': stop,
                'sequence': sequence,
            }


def bed_to_sequence_generator(
        bed_file: str,
        fasta_file: str,
        seq_as_str: bool = False):
    '''
    yields raw sequence info form a bed file and a fasta file.

    Args:
        bed_file: path to bed file.
        fasta_file: path to fasta file.
        seq_as_str: whether to set the 'sequence' value to be
            a python str type, or leave it as the Fasta Sequence type.
    '''

    bed_handle = BedTool(bed_file)

    fasta_handle = Fasta(fasta_file)

    for interval in bed_handle:
        chrom = interval.chrom
        start = interval.start
        stop = interval.stop

        sequence = fasta_handle[chrom][start:stop]
        if seq_as_str:
            sequence = str(sequence)
        yield {
            'chrom': chrom,
            'start': start,
            'stop': stop,
            'sequence': sequence,
        }





