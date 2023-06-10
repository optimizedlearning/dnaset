class BedLine():
	def __init__(self,s):
		s = s.split('\t')
		self.chrom = s[0]
		self.start = int(s[1])
		self.stop = int(s[2])
		
class RBedTool():
	'''
	A simple class for iterating over a .bed file
	'''
	def __init__(self,pth):
		#this is relatively fast, according to stack overflow
		with open(pth, "rbU") as f:
			self.sz = sum(1 for _ in f)
		
		self.fp = open(pth)
		self.line = 0

	def __getitem__(self,idx):
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
		self.line = idx+1
		return BedLine(fp.readline())

	def __iter__(self):
		return RBedToolIterator(self)

	def __len__(self):
		return self.sz

class RBedToolIterator:
	'''
	This is added so that multiple RBedToolIterators of the same RBedTool will not interfere with each other
	'''
	def __init__(self,rbedtool):
		self.rbedtool = rbedtool
		self.pos = 0
		self.it = 0

	def __next__(self):

		if self.it == len(self.rbedtool):
			raise StopIteration
		self.rbedtool.fp.seek(self.pos)
		self.rbedtool.line = self.it

		ret = self.rbedtool[self.it]

		self.pos = self.rbedtool.fp.tell()
		self.it+=1
		return ret

    

	
