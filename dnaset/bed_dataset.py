
import torch
from torch.utils.data import Dataset, IterableDataset
import pybedtools

class BedDataset(Dataset):
    '''
    very thin wrapper around a pybedtools.BedTool.
    '''
    def __init__(self, bed_file):
        if isinstance(bed_file, pybedtools.BedTool):
            self.bed =  bed_file
        else:
            self.bed = pybedtools.BedTool(bed_file)

    def __len__(self):
        return len(self.bed)


    def __getitem__(self, idx):
        return self.bed[idx]



