
import torch
from torch.utils.data import IterableDataset
import pyBigWig
from pyfaidx import Fasta
from pybedtools import BedTool
import tempfile
import random

from typing import Union, List, Optional, Callable

from dnaset import util
from dnaset.rbedtool import RBedTool


CATCH_GENERATOR_INDEX_ERROR = False

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
    enforce_common_length: bool = True,
    shuffle: Optional[bool] = False,
    chrom_ignore_chars: str = ""
        ):
    '''
    creates a bed file that tiles the provided genome, skipping gaps.
    each line of the output will be a sequence of length sequence_length.

    arguments:
        sequence_length: the sequence length for each line in the output.
        reference_fasta: the fasta file containing reference genome.
        gap_bed_list: a list of file names for the bed files containing
            assembly gaps.
        stride: distance between start of one sequence and start of next
            sequence (defaults to sequence_length).
        out_path: if specified, save the bed output here.
        enforce_common_length: if True, all output sequences must have length
            sequence_length. If False, there may be a few outputs with
            shorter lengths (because the areas between gaps in the gap files
            is not evenly divisible by sequence_length).
        shuffle: whether to shuffle the lines in the output bed file.
        chrom_ignore_chars: if specified, will exclude any chromosomes that include 
            any of the characters in the string 

    returns:
        a BedTool object for the resulting bed file.
        Note that if out_path is not supplied, the bed file is a temporary
            file.
        Probably the system has a policy to clean these up eventually, but it
            would be better for the  caller to take care of deleting it when
            finished otherwise.
    '''
    if stride is None:
        stride = sequence_length

    if isinstance(reference_fasta, str):
        reference_fasta = Fasta(reference_fasta)

    gap_bed_list = convert_single_to_list(gap_bed_list, str)
    if len(gap_bed_list) == 0:
        gap_bed = []
    else:
        gap_bed = BedTool(gap_bed_list[0])
        if len(gap_bed_list) > 1:
            gap_bed = gap_bed[0].cat(
                *[BedTool(bed)for bed in gap_bed_list[1:]])
    gap_bed = gap_bed.sort()

    fd, temp_file_name = tempfile.mkstemp(suffix='.bigwig_torch.bed')
    bed_fp = open(temp_file_name, mode='w+')

    assembly_starts = {
            chrom: 0 for chrom in reference_fasta.keys()
            }

    def write_bed(chrom, start, end):
        if any([chrom.find(c)!=-1 for c in chrom_ignore_chars]):
            return
        if start >= end:
            return
        if enforce_common_length and start + sequence_length > end:
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
    lines = []
    if shuffle:
        lines = bed_fp.readlines()
        random.shuffle(lines)
        bed_fp.seek(0)
        bed_fp.writelines(lines)



    # Open via file name so that it is inherited properly by child processes.
    out = BedTool(temp_file_name)
    if out_path is not None:
        out = out.moveto(out_path)
        bed_fp.close()

    # Note that since we made the tempfile with mkstemp, it will be the
    # caller's responsibility to cleanup the tempfile when out_path == None.
    # It may be ok most of the time to not clean it up and rely on some
    # system policy to eventually delete temporary files.

    return out




def bigwig_dataset_generator(
    bigwig_files: Union[List[str], str],
    reference_fasta: Union[str, Fasta],
    sequence_bed: Union[str, BedTool],
    sequence_transform: Callable = util.seq_to_array,
    start: int = 0,
    stop: int = -1,
    out_type: str = 'dict'
    ):
    '''
    generates numpy sequence arrays from a bigwig and bed file.

    arguments:
        bigwig_files: a single string or a list of strings specifying the paths to the bigwig files.
        reference_fasta: fasta file for the reference genome to use with these bigwig files.
        sequence_bed: a bed file specifying the intervals to use in the dataset.
        sequence_transform: a function called to transform the output sequence.
            defaults to converting the sequence into a numpy array.
        start: which row of the bed file to start at.
            Note that it will take O(start) time to yield the first value!
        stop: which row of the bed file to  end at (negative numbers index from the end of the file, just
            like regular array slicing)
	out_type: the output type. Can be 'dict' or 'tuple'
    '''

    bigwig_files = convert_single_to_list(bigwig_files, str)

    bigwigs = [pyBigWig.open(bw) for bw in bigwig_files]
    if not isinstance(sequence_bed, RBedTool):
        sequence_bed = RBedTool(sequence_bed)
    if not isinstance(reference_fasta, Fasta):
        reference_fasta = Fasta(reference_fasta, sequence_always_upper=True)
    
    def iterable_to_generator(iter):
        def gen(start,stop):
            if stop is None or stop==-1:
                stop = len(iter)
            while start < stop:
                yield iter[start]
                start += 1
        return gen

    def index_error_handler(msg):
        global CATCH_GENERATOR_INDEX_ERROR
        if CATCH_GENERATOR_INDEX_ERROR:
            print(msg)
            return None
        raise IndexError(msg)

    for interval in sequence_bed.slice(start,stop):
        chrom = interval.chrom
        start = interval.start
        stop = interval.stop
        # we copy to make the resulting array mutable

        if chrom is None: #handles parse error in rbedtool
            continue

        try:
            cause = "Fasta"
            sequence = sequence_transform(reference_fasta[chrom][start:stop])
            cause = "Bigwig"
            values = [bw.values(chrom, start, stop, numpy=True) for bw in bigwigs]
        except (IndexError):
            index_error_handler(
                f"Base pairs {start} to {stop} in chrom {chrom} are not accessible in {cause} file"
            )
            continue
        if out_type == 'dict':
            yield {
                    'sequence': sequence,
                    'values': values
                    }
        elif out_type == 'tuple':
            yield sequence, values
        else:
            raise ValueError(f"out type must be 'tuple' or 'int' (got {out_type})")

class catch_generator_index_error():
    '''
    Should only be used for with catch_out_of_range_generator():

    Will cause index errors raised by the main for-loop of bigwig_dataset_generator to be caught
    and printed to stderr
    '''
    def __enter__(self):
        global CATCH_GENERATOR_INDEX_ERROR
        CATCH_GENERATOR_INDEX_ERROR=True
    def __exit__(self, *args):
        global CATCH_GENERATOR_INDEX_ERROR
        CATCH_GENERATOR_INDEX_ERROR=True



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
        input_bed_file: str,
        sequence_transform: Callable= util.seq_to_array,
        out_type = "dict"
    ):
        '''
        arguments:
            bigwig_files: a single string or a list of strings specifying the paths to the bigwig files.
            reference_fasta_file: fasta file for the reference genome to use with these bigwig files.
            input_bed_file: a bed file path specifying the intervals to use in the dataset
            sequence_transform: a function called to transform the output sequence.
                defaults to converting the sequence into a numpy array.
            out_type: format of the output.
                dict -> __getitem__ will return a dict with fields 'sequence' and 'values'
                tuple -> __getitem__ will return a tuple (sequence, values)
                other -> will raise error
        '''

        self.reference_fasta_file = reference_fasta_file
        self.reference_fasta = Fasta(
            reference_fasta_file,
            sequence_always_upper=True)

        self.input_bed_file = input_bed_file
        self.input_bed = RBedTool(input_bed_file)

        self.sequence_transform = sequence_transform

        # if we need random access to the bed file, this attribute will
        # be populated because BedTools[idx] is O(idx).
        self.random_access_input_bed = None

        self.length = None

        bigwig_files = convert_single_to_list(bigwig_files, str)

        self.bigwig_files = bigwig_files

        self.open_bigwig_files = [pyBigWig.open(bw) for bw in bigwig_files]

        supported_out_types = ["dict","tuple"]
        if out_type not in supported_out_types:
            raise ValueError(f"out type must be 'tuple' or 'dict' (got {out_type})")
        self.out_type = out_type

    def __len__(self):
        # rbedtools will compute the length of the file in lines at initialization 
        return len(self.input_bed)

    def __getitem__(self, idx):
        

        bed_line = self.input_bed[idx]

        chrom = bed_line.chrom
        start = bed_line.start
        stop = bed_line.stop

        # we copy to make the resulting array mutable
        sequence = self.sequence_transform(
            self.reference_fasta[chrom][start:stop]
        )

        values = [
            bw.values(chrom, start, stop, numpy=True)
            for bw in self.open_bigwig_files
        ]
        if self.out_type=='dict':
            return {
                    'sequence': sequence,
                    'values': values
                    }
        else:
            return sequence, values

    def __iter__(self):
        worker_info = torch.utils.data.get_worker_info()
        if worker_info is None:
            iter_start = 0
            iter_end = self.length
        else:
            # following line uses the identity ceil x/y = floor (y+x-1)/y
            num_workers = worker_info.num_workers
            lines_per_worker = (self.length + num_workers - 1) // num_workers

            iter_start = worker_info.worker_id * lines_per_worker
            iter_end = min(iter_start + lines_per_worker, self.length)

        return bigwig_dataset_generator(
                    self.bigwig_files,
                    self.reference_fasta_file,
                    self.input_bed_file,
                    self.sequence_transform,
                    iter_start,
                    iter_end,
                    self.out_type
                )

