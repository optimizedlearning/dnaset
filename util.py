import numpy as np


# class SeqArrayWrapper:
#     def __init__(self, s):
#         self.s = memoryview(str(s).encode('ascii'))
#     @property
#     def __array_interface__(self):
#         return {
#             'shape': (len(self.s),),
#             'typestr': '|u1',
#             'version': 3,
#             'data': self.s
#         } 

# def seq_to_array(s):
#     a = np.asarray(SeqArrayWrapper(s))
#     a.setflags(write=True)
#     return a

def seq_to_array(s):
    return np.asarray(memoryview(str(s).encode('ascii')), dtype=np.uint8)


def char_view(s):
    return s.view(dtype='|S1')

def int_view(s):
    return s.view(dtype=np.uint8)
