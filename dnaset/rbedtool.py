CATCH_BED_PARSE_ERROR = False
current_line = 0

class BedFrame():
	def __init__(self,s):
		global current_line
		global CATCH_BED_PARSE_ERROR
		try:
			a = s.split('\t')
			self.chrom = a[0]
			self.start = int(a[1])
			self.stop = int(a[2])
		except:
			s = f"There was an issue while trying to parse BedFrame input\nline {current_line}: {s}"
			if CATCH_BED_PARSE_ERROR:
				print(s)
				print("catch_bed_parser_error is enabled. self.chrom will be None")
				self.chrom = None
				self.start = None
				self.stop = None
			else:
				raise ValueError(s)

class RBedTool():
	'''
	A simple class for iterating over the lines of a bed file.
	Accessing elements sequentially (such as with an iterator) is
	O(1). Otherwise, the BedTool will perform a linear search
	for the accessed line.

	__getitem__() will return a BedFrame object, a simple struct with
	feilds chrom, stop, and start
	'''
	def __init__(self,pth):
		'''
			args:
				pth - the path to your bed file
		'''
		#this is relatively fast, according to stack overflow
		with open(pth, "rbU") as f:
			self.sz = sum(1 for _ in f)
		
		self.fp = open(pth)
		self.line = 0
		self.sl = None

	def __getitem__(self,idx):
		global current_line

		if type(idx)==slice:
			raise TypeError("Slice passed to RBedTool.__getitem__, please use RBedTool.slice for slicing")
		def line_align(fp):
			#backs up a file pointer character-by-character until it reaches a newline
			while fp.tell()!=0 and fp.read(1)!='\n':
				fp.seek(fp.tell()-2)
		
		if idx >= len(self):
			raise IndexError
		j = 64
		fp = self.fp
		while self.line > idx and self.line!=0:
			while fp.tell() < j-1 and j > 1:
				j = j//2
			fp.seek(fp.tell()-j)
			self.line -= fp.read(j).count('\n')
			fp.seek(fp.tell()-j)
		line_align(fp)

		for _ in range(idx-self.line):
			fp.readline()

		current_line = self.line
		self.line = idx+1

		return BedFrame(fp.readline())

	def __iter__(self):
		return RBedToolIterator(self)

	def __len__(self):
		return self.sz

	def slice(self,start = None, stop = None, step = None):
		'''
		Returns an iterator for this RBedTool over a slice of the data.
		As far as I can tell, itertools does not support lazy loading for 
		sliced data, so I implemented a separate method
		args:
			sl - a slice object (e.g. 10:24:2, 1:5, etc.)
		returns:
			RBedTool iterator that will iterate over the values specified by the slice
		'''
		sl = {
			'start': start,
			'stop': stop,
			'step': step
		}
		if sl['start'] is None:
			sl['start'] = 0
		if sl['stop'] is None:
			sl['stop'] = len(self)
		if sl['step'] is None:
			sl['step'] = 1

		self.sl = sl
		return self

class RBedToolIterator:
	'''
	This is added so that multiple RBedToolIterators of the same RBedTool will not interfere with each other
	'''
	def __init__(self,rbedtool,sl=None):
		self.rbedtool = rbedtool
		sl = rbedtool.sl



		if sl is not None:
			self.line = sl['start'] 
			self.stop = len(rbedtool) if sl['stop']==-1 else sl['stop'] #the line this iterator will stop at
			self.step = sl['step']
		else:
			self.line = 0
			self.stop = len(rbedtool)
			self.step = 1
		self.pos = 0
		self.init = False #Will be used to determine if rbedtool.line should be overridden on a given iteration
		
		#to prevent future calls to rbedtool.__iter__() from usuing this slice
		rbedtool.sl = None

	def __next__(self):
		fp = self.rbedtool.fp
		rbedtool = self.rbedtool

		if self.line == self.stop:
			raise StopIteration
		
		if self.init: #We don't want to call fp.seek() until we know where line self.line is
			fp.seek(self.pos)
			rbedtool.line = self.line
			for _ in range(self.step):
				if self.line==self.stop:
					raise StopIteration
				ret = rbedtool[self.line]
				self.line+=1
		else:
			ret = rbedtool[self.line]
			self.line+=1
		
		self.pos = fp.tell()
		self.init = True

		return ret

    

	
class catch_bed_parse_error():
	def __enter__(self):
		global CATCH_BED_PARSE_ERROR
		CATCH_BED_PARSE_ERROR = True
	def __exit__(self,*args):
		global CATCH_BED_PARSE_ERROR
		CATCH_BED_PARSE_ERROR = False