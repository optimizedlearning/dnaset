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
    array = torch.Tensor(memoryview(str(s).encode('ascii'))) #I'm not entirely clear on what memoryview is used for
    array = array.type(torch.uint8)
    if copy:
        array = torch.clone(array)
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

'''
def bed_to_bbed(pth, out_pth, n_chroms = None):

    Converts a bed file into a compressed format that allows for random line access

    args:
        pth - location of bed file
        n_chroms - number of chromosomes expected in the BED file. Will throw an error
            if unique chromosomes exceeds this number

    if n_chroms is not None:
        chroms = [None for _ in range(n_chroms)]
    else:
        chroms = []
    i = 0
    with open(pth) as ifp, open(out_pth,'wb') as ofp:
        for line in ifp:
            chrom = line.split('\t')[0]
            if chrom not in chroms:
                if n_chroms is None:
                    chroms.append(chrom)
                else:
                    try:
                        chroms[i]=chrom
                    except(IndexError):
                        raise AssertionError(f"Number of unique chroms exceeds expected value ({n_chroms})\nChromosomes identified: {chroms}")

    
'''





