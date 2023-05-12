
import torch
from torch.utils.data import Dataset, IterableDataset
import pyBigWig
from pybedtools import BedTool
from pyfaidx import Fasta
import numpy as np
import tempfile
import os
import itertools

from typing import Any, Union, List, Optional

import util

# hack to memorize an iterator so that we can do some faster random access.
# obviously only works if you have enough memory to hold everything.
# TODO: do some fancy ness with mutiple iterators to make this better.
# Right now, this is not really any better than list(iter(base))
class RandomAccessViaIter:

    def __init__(self, base):
        self.base = base
        self.it = enumerate(self.base)
        self.max_index = -1
        self.cached_values = {
        }

    def __getitem__(self, idx):
        if idx > len(base) or idx < 0:
            raise IndexError
        if idx not in cached_values:
            while self.max_index < idx:
                self.max_index, value = next(self.it)
                self.cached_values[self.max_index] = value

        return cached_values[idx]

    

def tile_genome(
    sequence_length: int,
    reference_fasta: Union[str, Fasta],
    gap_bed_list: Union[List[Union[BedTool, str]], BedTool, str],
    stride: Optional[int] = None,
    out_path: Optional[str] = None,
    shuffle: Optional[bool] =False):
    '''
    creates a bed file that tiles the provided genome, skipping gaps.
    each line of the output will be a sequence of length sequence_length.

    arguments:
        sequence_length: the sequence length for each line in the output.
        reference_fasta: the fasta file containing reference genome.
        gap_bed_list: a list of file names for the bed files containing assembly gaps.
        stride: distance between start of one sequence and start of next sequence (defauots to sequence_length).
        out_path: if specified, save the bed output here.
        shuffle: whether to shuffle the lines in the output bed file.

    returns:
        a BedTool object for the resulting bed file.
        Note that if out_path is not supplied, the bed file is a  temporary file.
        Probably the system has a policy to clean these up eventually, but it would be
        better for the  caller to take care of deleting it when finished otherwise.
    '''
    if stride is None:
        stride = sequence_length

    if not isinstance(reference_fasta, Fasta):
        reference_fasta = Fasta(reference_fasta)
    if isinstance(gap_bed_list, str):
        gap_bed_list = [gap_bed_list]


    if len(gap_bed_list) == 0:
        gap_bed = []
    else:
        gap_bed  = BedTool(gap_bed_list[0])
        if len(gap_bed_list) > 1:
            gap_bed = gap_bed[0].cat(*[BedTool(bed)for bed in gap_bed_list[1:]])
    gap_bed = gap_bed.sort()

    fd, temp_file_name = tempfile.mkstemp(suffix='.bigwig_torch.bed')
    bed_fp = open(temp_file_name, mode='w')

    assembly_starts = {
            chrom: 0 for chrom in reference_fasta.keys()
            }

    def write_bed(chrom, start, end):
        if start >= end:
            return

        for line_start in range(start, end-start, stride):
            bed_fp.write(
                    f"{chrom}\t{line_start}\t{line_start + sequence_length}\n"
                    )


    for gap in gap_bed:
        chrom = gap.chrom
        gap_start = gap.start
        gap_end = gap.end

        write_bed(chrom, assembly_starts[chrom], gap_start)

        assembly_starts[chrom] = gap_end

    for chrom in reference_fasta.keys():
        write_bed(chrom, assembly_starts[chrom], len(reference_fasta[chrom]))

    bed_fp.seek(0)

    if shuffle:
        lines = bed_fp.readlines()
        random.shuffle(lines)
        bed_fp.seek(0)
        bed_fp.write(lines)
        bed_fp.seek(0)

    # Open via file name so that it is inhereted properly by child processes.
    out = BedTool(temp_file_name)
    if out_path is not None:
        bed_file_name  = out_path
        out = out.moveto(out_path)
        bed_fp.close()

    return out

def bigwig_dataset_generator(
    bigwig_filenames,
    reference_fasta,
    sequence_bed,
    start=0,
    end=-1):
    '''
    generates numpy sequence arrays from a bigwig and bed file.
    
    arguments:
        bigwig_files: a single string or a list of strings specifying the paths to the bigwig files.
        reference_fasta: fasta file for the reference genome to use with these bigwig files.
        sequence_bed: a bed file specifying the intervals to use in the dataset.
        start: which row of the bed file to start at.
            Note that it will take O(start) time to yield the first value!
        end: which row of the bed file to  end at (negative numbers index from the end of the file, just
            like regular array slicing)
    '''

    if isinstance(bigwig_filenames, str):
        bigwig_filenames = [bigwig_filenames]
    bigwigs = [pyBigWig.open(bw) for bw in bigwig_filenames]
    if not isinstance(sequence_bed, BedTool):
        sequence_bed = BedTool(sequence_bed)
    if not isinstance(reference_fasta, Fasta):
        reference_fasta = Fasta(reference_fasta, sequence_always_upper=True)

    for interval in sequence_bed[start:stop]:
        chrom = bed_line.chrom
        start = bed_line.start
        stop = bed_line.stop

        sequence = np.asarray(reference_fasta[chrom][start:stop])
        values = [bw.values(chrom, start, stop, numpy=True) for bw in bigwigs]

        yield {
                'sequence': sequence,
                'values': values
                }

class BigWigDataset(Dataset):
    '''
    a wrapper for porting bigwig tracks to a pytorch dataset.

    the outputs are presented as dictionaries:
    {
        'sequence': a numpy array of characters representing the genome sequence
        'values': a list of numpy arrays containing the corresponding values.
    }
    '''
    def __init__(self, bigwig_files: Union[List[str], str], reference_fasta, sequence_bed=None, sequence_length=None, gap_bed=None, stride=None, out_path=None):
        '''
        arguments:
            bigwig_files: a single string or a list of strings specifying the paths to the bigwig files.
            reference_fasta: fasta file for the reference genome to use with these bigwig files.
            sequence_bed: a bed file specifying the intervals to use in the dataset (optional; if not specified, then
                sequence_length must be specified).
            sequence_length: the length of the output sequences to extract from the reference_fasta. Should only be
                used when sequence_bed is None. The dataset will tile the reference genome with sequences of length sequence_length.
            gap_bed: if sequence_length is not None, this is a bed file specifying gaps that  will not be tiled in the reference.
            stride: if sequence_length is not None, this specifies a stride for tiling the reference. If unspecified, will default to
                sequence_length.
            out_path: if sequence_bed is not specified, save the generated sequence_bed here (no saving if None)
        '''

        self.do_delete_bed  = out_path is None and sequence_bed is None

        if sequence_length is not None and sequence_bed is not None:
            raise ValueError('Must specify exactly one of sequence_length and sequence_bed')
        if sequence_length is None and sequence_bed is None:
            raise ValueError('Must specify exactly one of sequence_length and sequence_bed')

        self.reference_fasta = Fasta(reference_fasta, sequence_always_upper=True)

        self.sequence_length = sequence_length
        if stride is None:
            stride  = sequence_length
        self.stride =  stride

        

        if sequence_bed is not None:
            self.sequence_bed = BedTool(sequence_bed)
        elif sequence_length is not None:
            if isinstance(gap_bed, str):
                gap_bed = [gap_bed]
            self.sequence_bed = tile_genome(sequence_length, reference_fasta, gap_bed, stride, out_path, shuffle=True)

        self.random_access_sequence_bed = None

        self.length = len(self.sequence_bed)



        if isinstance(bigwig_files, str):
            bigwig_files = [bigwig_files]

        self.bigwig_files = bigwig_files

        self.open_bigwig_files = [pyBigWig.open(bw) for bw in bigwig_files]



    def __len__(self):
        return self.length

    def __del__(self):

        if self.do_delete_bed:
            # we do this in a try/catch because the file might be deleted by another process
            # if we are doing parallel data loaders.
            try:
                os.unlink(self.sequence_bed.fn)
            except FileNotFoundError:
                pass


    def __getitem__(self, idx):

        # bed_file[idx] is actually an O(idx) time operation.
        # so, we need some kind of hack to speed things up.
        if self.random_access_sequence_bed is None:
            self.random_access_sequence_bed = list(iter(self.sequence_bed))

        bed_line = self.random_access_sequence_bed[idx]

        chrom = bed_line.chrom
        start = bed_line.start
        stop = bed_line.stop

        # we copy to make the resulting array mutable
        sequence = util.seq_to_array(self.reference_fasta[chrom][start:stop]).copy()
        values = [bw.values(chrom, start, stop, numpy=True) for bw in self.open_bigwig_files]

        return {
                'sequence': sequence,
                'values': values
                }


    def __iter__(self):
        worker_info = torch.utils.data.get_worker_info()
        if worker_info is None:
            iter_start = 0
            iter_end = self.length
        else:
            # floor (y+x-1)/y = ceil x/y
            lines_per_worker = (self.length + worker_info.num_workers - 1) // worker_info.num_workers

            iter_start = worker_info.worker_id * lines_per_worker
            iter_end = min(iter_start + lines_per_worker, self.length)

        return iter(
            bigwig_dataset_generator(
                self.bigwig_files,
                self.reference_fasta,
                self.sequence_bed,
                iter_start, iter_end
            )
        )




