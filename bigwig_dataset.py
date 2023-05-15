
import torch
from torch.utils.data import Dataset, IterableDataset
import pyBigWig
from pybedtools import BedTool
from pyfaidx import Fasta
import numpy as np
import tempfile

from typing import Any, Union, List, Optional

import util


def convert_single_to_list(maybe_list, base_type):
    if isinstance(maybe_list, base_type):
        return [maybe_list]
    return maybe_list


def tile_genome(
    sequence_length: int,
    reference_fasta: Union[str, Fasta],
    gap_bed_list: Union[List[str], str],
    stride: Optional[int] = None,
    out_path: Optional[str] = None,
    shuffle: Optional[bool] = False):
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
        Note that if out_path is not supplied, the bed file is a temporary file.
        Probably the system has a policy to clean these up eventually, but it would be
        better for the  caller to take care of deleting it when finished otherwise.
    '''
    if stride is None:
        stride = sequence_length

    if isinstance(reference_fasta, str):
        reference_fasta = Fasta(reference_fasta)

    gap_bed_list = convert_single_to_list(gap_bed_list, str)
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

    # Open via file name so that it is inherited properly by child processes.
    out = BedTool(temp_file_name)
    if out_path is not None:
        bed_file_name  = out_path
        out = out.moveto(out_path)
        bed_fp.close()

    # Note that since we made the tempfile with mkstemp, it will be the
    # caller's responsibility to cleanup the tempfile when out_path == None.
    # It may be ok most of the time to not clean it up and rely on some system policy.

    return out

def bigwig_dataset_generator(
    bigwig_files: Union[List[str], str],
    reference_fasta: Union[str, Fasta],
    sequence_bed: Union[str, BedTool],
    start=0,
    stop=-1):
    '''
    generates numpy sequence arrays from a bigwig and bed file.

    arguments:
        bigwig_files: a single string or a list of strings specifying the paths to the bigwig files.
        reference_fasta: fasta file for the reference genome to use with these bigwig files.
        sequence_bed: a bed file specifying the intervals to use in the dataset.
        start: which row of the bed file to start at.
            Note that it will take O(start) time to yield the first value!
        stop: which row of the bed file to  end at (negative numbers index from the end of the file, just
            like regular array slicing)
    '''

    bigwig_files = convert_single_to_list(bigwig_files, str)

    bigwigs = [pyBigWig.open(bw) for bw in bigwig_files]
    if not isinstance(sequence_bed, BedTool):
        sequence_bed = BedTool(sequence_bed)
    if not isinstance(reference_fasta, Fasta):
        reference_fasta = Fasta(reference_fasta, sequence_always_upper=True)
    if not isinstance(sequence_bed, BedTool):
        sequence_bed = BedTool(sequence_bed)

    for interval in sequence_bed[start:stop]:
        chrom = interval.chrom
        start = interval.start
        stop = interval.stop

        # we copy to make the resulting array mutable
        sequence = util.seq_to_array(reference_fasta[chrom][start:stop]).copy()
        values = [bw.values(chrom, start, stop, numpy=True) for bw in bigwigs]

        yield {
                'sequence': sequence,
                'values': values
                }

class BigWigDataset(IterableDataset):
    '''
    a wrapper for porting bigwig tracks to a pytorch dataset.

    the outputs are presented as dictionaries:
    {
        'sequence': a numpy array of ascii codes representing the genome sequence
        'values': a list of numpy arrays containing the corresponding values.
    }
    '''
    def __init__(
        self,
        bigwig_files: Union[List[str], str],
        reference_fasta_file: str,
        input_bed_file: str
    ):
        '''
        arguments:
            bigwig_files: a single string or a list of strings specifying the paths to the bigwig files.
            reference_fasta_file: fasta file for the reference genome to use with these bigwig files.
            input_bed_file: a bed file path specifying the intervals to use in the dataset
        '''

        self.reference_fasta_file = reference_fasta_file
        self.reference_fasta = Fasta(reference_fasta_file, sequence_always_upper=True)


        self.input_bed_file = input_bed_file
        self.input_bed = BedTool(input_bed_file)

        # if we need random access to the bed file, this attribute will
        # be populated because BedTools[idx] is O(idx).
        self.random_access_input_bed = None

        self.length = len(self.input_bed)


        bigwig_files = convert_single_to_list(bigwig_files, str)

        self.bigwig_files = bigwig_files

        self.open_bigwig_files = [pyBigWig.open(bw) for bw in bigwig_files]



    def __len__(self):
        return self.length



    def __getitem__(self, idx):
        # bed_file[idx] is actually an O(idx) time operation.
        # so, we need some kind of hack to speed things up.
        if self.random_access_input_bed is None:
            self.random_access_input_bed = list(iter(self.input_bed))

        bed_line = self.random_access_input_bed[idx]

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
        print("ran iter")
        worker_info = torch.utils.data.get_worker_info()
        if worker_info is None:
            iter_start = 0
            iter_end = self.length
        else:
            # following line uses the identity ceil x/y = floor (y+x-1)/y
            lines_per_worker = (self.length + worker_info.num_workers - 1) // worker_info.num_workers

            iter_start = worker_info.worker_id * lines_per_worker
            iter_end = min(iter_start + lines_per_worker, self.length)

        return bigwig_dataset_generator(
                   self.bigwig_files,
                   self.reference_fasta_file,
                   self.input_bed_file,
                   iter_start, iter_end
                )




