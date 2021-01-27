"""
algs
=================================================
The core algorithm module of my BMI 203 Project 1
"""

import numpy as np

class PairwiseAligner:
	"""Paiwise Alignment Base Class

	This is the base (parent) class for all pairwise alignment algorithms
	implemented in this API. All pariwise alignment algorithm classes
	inheret from this base class. This class can be used to create
	any pairwise alignment algorithm class.

	Parameters
	----------
	substitution_matrix_file : str
		path/filename of substitution matrix
    gap_open : float
		gap opening penalty (should be positive value, think abs value)
    gap_extend : float
	 	gap extension penalty (should be postivie value, think abs value)

	Attributes
	----------
	seq1 : str
		sequence from fasta1
	seq2 : str
		sequence from fasta2
	header1 : str
		fasta header from fasta1
	header2 : str
		fasta header from fasta2
	alignment_ : tuple
		Tuple of format (seq1_alignment, seq2_alignment)
	alignment_score_ : float
		Score of alignment from algorithm

	"""
	def __init__(self, substitution_matrix_file, gap_open, gap_extend):
		# generating substitution matrix
		self.sub_matrix = self._initiate_sub_matrix(substitution_matrix_file)
		self.gap_open = gap_open
		self.gap_extend = gap_extend
		self.seq1 = None
		self.seq2 = None
		self.score_matrix = None

	def _initiate_sub_matrix(self, matrix_file):
		"""This function reads in a scoring matrix from any matrix like file.
		Where there is a line of the residues followed by substitution
		matrix.

		Parameters
		----------
		matrix_file : str
			name (and associated path if not in current
			working directory) of the matrix file that contains the
			scoring matrix.

		Returns
		-------
		dict_sm : dict
		 	substitution matrix dictionary with tuple
			of the two residues as the key and score as value
			e.g. {('A', 'A') -> 4} or {('A', 'D') -> -8}

		"""
		with open(matrix_file, 'r') as f:
			dict_sm = {} # Dictionary for storing scores
			residue_list = [] # For storing residue list
			start = False # trigger for reading in score values
			res_2 = 0 # used for generating substitution matrix
			# reading file line by line
			for line_num, line in enumerate(f):
	            # Reading in Residue List
				if '#' not in line.strip() and start == False:
					residue_list = [k for k in line.strip().upper().split(' ') if k != '']
					start = True
				# Generating substitution scoring dictionary
				elif start == True and res_2 < len(residue_list):
					line = [k for k in line.strip().split(' ') if k != '']
					# line.strip().replace('  ', ' ').split(' ')
					# reading in line by line to create substitution dictionary
					assert len(residue_list) == len(line)
					for res_1 in range(len(line)):
						dict_sm[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
						# dict_sm[(residue_list[res_2], residue_list[res_1])] = float(line[res_1])
					res_2 += 1
				elif start == True and res_2 == len(residue_list):
					break
			return dict_sm

	def _read_fasta(self, fasta_file):
		"""This function reads in a FASTA file and returns the associated
		string of characters (residues or nucleotides) and the header.

		Parameters
		----------
		fasta_file : str
			name (and associated path if not in current working directory)
			of the Fasta file.

		Returns
		-------
		seq : str
			String of characters from FASTA file
		header : str
			Fasta header

		"""
		with open(fasta_file) as f:
			is_header = True
			seq = '' # initializing sequence
			for line in f:
				if is_header and line.strip().startswith('>'):
					header = line.strip() # reading in fasta header
					is_header = False
				else:
					seq += line.strip().upper() # generating full seq
		return seq, header

	def _initiate_seqs(self, fasta_1, fasta_2):
		"""This function initiates the sequences for alignment by
		calling the _read_fasta method on the two fasta files.
		This saves the sequences and headers into the fields
		self.seq1, self.seq2, self.header1, self.header2.

		Parameters
		----------
		fasta_1 : str
			path and filename of first fasta for alignment.
		fasta_2 : str
			path and filename of second fasta for alignment.

		"""
		seq1, header1 = self._read_fasta(fasta_1)
		seq2, header2 = self._read_fasta(fasta_2)
		self.seq1 = seq1
		self.seq2 = seq2
		self.header1 = header1
		self.header2 = header2

	def _initiate_matrix(self):
		"""This function initializes the score matrix of dimensions
		len(seq1)+1 x len(seq2)+1. This method must be run after the
		_initate_seqs method.

		Returns
		-------
		score_matrix : array-like
			array of size len(seq1)+1 x len(seq2)+1

		"""
		seq1_len = len(self.seq1) + 1
		seq2_len = len(self.seq2) + 1
		score_matrix = np.zeros((seq1_len, seq2_len))
		return score_matrix

	def _initate_all_matrices(self):
		"""This function initates all of the score matrices to the formats
		that are not specific to a given alignment algorithm.

		Returns
		-------
		s_mat : array-like
			initialized core scoring matrix
		x_mat : array-like
			initialized X gap scoring matrix, seq2 gaps
		y_mat : array-like
		 	initialized Y gap scoring matrix, seq1 gaps
		tb_mat : array-like
			initalized traceback matrix

		"""
		assert self.seq1 != None
		assert self.seq2 != None
		seq1 = self.seq1
		seq2 = self.seq2
		# Initializing Core Score Matrix
		self.score_matrix = self._initiate_matrix()
		s_mat = self.score_matrix
		# Initializing Additional Scoring Matrices for Indels (Gaps)
		x_mat = self._initiate_matrix() # seq2 gap matrix
		x_mat[:, 0] = np.arange(len(seq1)+1)*-self.gap_extend - self.gap_open
		x_mat[0, :] = np.arange(len(seq2)+1)*-self.gap_open*100000000
		y_mat = self._initiate_matrix() # seq1 gap matrix
		y_mat[:, 0] = np.arange(len(seq1)+1)*-self.gap_open*100000000
		y_mat[0, :] = np.arange(len(seq2)+1)*-self.gap_extend - self.gap_open
		# Initialize Traceback matrix
		tb_mat = self._initiate_matrix().astype('str')
		return s_mat, x_mat, y_mat, tb_mat

	def _fill_scoring_matrices(self, s_mat, x_mat, y_mat, tb_mat):
		"""This function must be overwritten with a new method
		for filling in the scoring matrices using affine gap
		for either local or global alignment, depending on
		the alignment algorithm.

		:meta private:
		:noindex:

		"""
		pass

	def _traceback(self, s_mat, tb_mat):
		"""This function must be overwritten with a new method
		for tracing back through the scoring matrices for a
		given alignment algorithm.

		:meta private:
		:noindex:

		"""
		pass

	def _algorithm_matrix_initiation(self, s_mat):
		"""This function must be overwtitten with a method for
		initializing the core score matrix, or s_mat in this
		class implementation, such that it is initalized for
		the given alignment algorithm. This differs for local
		and global alignment.

		:meta private:
		:noindex:

		"""
		pass

	def align(self, fasta1, fasta2, return_alignment=True):
		"""This is a method that must be overwritten when inhereting from
		this class. for NeedlemanWunsch

		Parameters
		----------
		fasta1 : str
			Filename string of first fasta for alignment
		fasta2 : str
			Filename string of second fasta for alignment
		return_alignment : bool, default=True
			Whether or not to return the alignmenti if True then
			the alignment is returned and if False then None
			is returned

		Returns
		-------
		alignment : tuple
			Tuple of format (seq1_alignment, seq2_alignment)

		"""
		# Initializing Sequences
		self._initiate_seqs(fasta1, fasta2)
		#  Initializing All Matrices
		s_mat, x_mat, y_mat, tb_mat = self._initate_all_matrices()
		s_mat = self._algorithm_matrix_initiation(s_mat)
		s_mat, x_mat, y_mat, tb_mat = self._fill_scoring_matrices(s_mat, x_mat, y_mat, tb_mat)
		if return_alignment == True:
			alignment, alignment_score = self._traceback(s_mat, tb_mat)
			return alignment


################################################################################


class SmithWaterman(PairwiseAligner):
	"""This class is an implenetation of the Smith Waterman pairwise
	local aligner algorithm and inherets from the PairwiseAligner
	class. This class can be used to align any two protein sequences
	for SmithWaterman.

	Parameters
	----------
	substitution_matrix_file : str
		path/filename of substitution matrix
	gap_open : float
		gap opening penalty
	gap_extend : float
		gap extension penalty

	Attributes
	----------
	seq1 : str
		sequence from fasta1
	seq2 : str
		sequence from fasta2
	header1 : str
		fasta header from fasta1
	header2 : str
		fasta header from fasta2
	alignment_ : tuple
		Tuple of format (seq1_alignment, seq2_alignment)
	alignment_score_ : float
		Score of alignment from algorithm

	"""
	def __init__(self, substitution_matrix_file, gap_open, gap_extend):
		super(SmithWaterman, self).__init__(substitution_matrix_file, gap_open, gap_extend)

	def _fill_scoring_matrices(self, s_mat, x_mat, y_mat, tb_mat):
		"""This function uses dynamic programming to fill in all of the
		scoring and traceback matrices. This specific implementation
		is local alignment using the SmithWaterman alignment algorithm.
		The function takes in all initalized scoring matrices and fills
		them out for SmithWaterman.

		Parameters
		----------
		s_mat : array-like
			initialized core scoring matrix
		x_mat : array-like
			initialized X gap scoring matrix, seq2 gaps
		y_mat : array-like
			initialized Y gap scoring matrix, seq1 gaps
		tb_mat : array-like
			initalized traceback matrix

		Returns
		-------
		s_mat : array-like
			Filled core scoring matrix
		x_mat : array-like
			Filled X gap scoring matrix, seq2 gaps
		y_mat : array-like
			Filled Y gap scoring matrix, seq1 gaps
		tb_mat : array-like
			Filled Traceback matrix, used for constructing alignment

		"""
		# Filling in score Matrix
		seq1 = self.seq1
		seq2 = self.seq2
		for i in range(s_mat.shape[0])[1:]:
			for j in range(s_mat.shape[1])[1:]:
				# t gap matrix filling
				gap_dict = {1: 'new_gap', 0: 'extended_gap'}
				new_gap = (s_mat[i-1, j] - self.gap_open)
				extended_gap = (x_mat[i-1, j] - self.gap_extend)
				x_mat[i, j] = max([extended_gap, new_gap])
				gap_type_x = gap_dict[np.argmax([extended_gap, new_gap])]
				# u gap matrix filling
				new_gap = (s_mat[i, j-1] - self.gap_open)
				extended_gap = (y_mat[i, j-1] - self.gap_extend)
				y_mat[i, j] = max([extended_gap, new_gap])
				gap_type_y = gap_dict[np.argmax([extended_gap, new_gap])]
				# Taking maximum transition
				top = x_mat[i, j]
				side = y_mat[i, j]
				diag = s_mat[i-1, j-1] + self.sub_matrix[(seq1[i-1], seq2[j-1])]
				s_mat[i,j] = max((top, side, diag, 0))
				max_ind = np.argmax((top, side, diag, 0))
				# filling in traceback matrix
				if max_ind == 0 and gap_type_x == 'new_gap':
					tb_mat[i, j] = 'top_new'
				elif max_ind == 0 and gap_type_x == 'extended_gap':
					tb_mat[i, j] = 'top_extended'
				elif max_ind == 1 and gap_type_y == 'new_gap':
					tb_mat[i, j] = 'side_new'
				elif max_ind == 1 and gap_type_y == 'extended_gap':
					tb_mat[i, j] = 'side_extended'
				elif max_ind == 2:
					tb_mat[i, j] = 'diag'
				elif max_ind == 3:
					tb_mat[i, j] = 'end'
		self.align_score_ = np.amax(s_mat)
		return s_mat, x_mat, y_mat, tb_mat

	def _traceback(self, s_mat, tb_mat):
		"""This functions uses the core scoring matrix and the traceback
		matrix to find the alignment score and full alignment of seq1
		and seq2. This function is implemented for the local alignment
		SmithWaterman algorithm. for SmithWaterman

		Parameters
		----------
		s_mat : array-like
			Core scoring matrix
		tb_mat : array-like
			Traceback matrix, used for constructing alignment

		Returns
		-------
		alignment : tuple
			Tuple of format (seq1_alignment, seq2_alignment)
		alignment_score : float
			Float value of maximum alignment score in the core score matrix

		"""
		# Traceback
		seq1 = self.seq1
		seq2 = self.seq2
		alignment_seq1 = ''
		alignment_seq2 = ''
		max_ind = np.unravel_index(np.argmax(s_mat, axis=None), s_mat.shape)
		i, j = max_ind
		while i > 0 and j > 0:
			tb_val = tb_mat[i, j]
			if tb_val == 'top_extended':
				alignment_seq2 = '-' + alignment_seq2
				alignment_seq1 = seq1[i-1] + alignment_seq1
				i -= 1
			elif tb_val == 'top_new':
				alignment_seq2 = '-' + alignment_seq2
				alignment_seq1 = seq1[i-1] + alignment_seq1
				i -= 1
			elif tb_val == 'side_extended':
				alignment_seq2 = seq2[j-1] + alignment_seq2
				alignment_seq1 = '-' + alignment_seq1
				j -= 1
			elif tb_val == 'side_new':
				alignment_seq2 = seq2[j-1] + alignment_seq2
				alignment_seq1 = '-' + alignment_seq1
				j -= 1
			elif tb_val == 'diag':
				alignment_seq2 = seq2[j-1] + alignment_seq2
				alignment_seq1 = seq1[i-1] + alignment_seq1
				i -= 1
				j -= 1
			elif tb_val == 'end':
				break
		self.alignment_ = (alignment_seq1, alignment_seq2)
		alignment = self.alignment_
		alignment_score = self.align_score_
		return alignment, alignment_score


	def _algorithm_matrix_initiation(self, s_mat):
		"""This function intialized the core score matrix for the
		SmithWaterman Alignment by setting the first row and
		column to be filled with zeros, because the local
		alignment does not allow for negative values.
		for SmithWaterman

		Parameters
		----------
		s_mat : array-like
			Initialized core scoring matrix

		Returns
		-------
		s_mat : array-like
			Core scoring matrix initalized for SmithWaterman

		"""
		# Filling in first column and row
		s_mat[1:, 0] = 0
		s_mat[0, 1:] = 0
		return s_mat


class NeedlemanWunsch(PairwiseAligner):
	"""This class is an implenetation of the NeedlemanWunsch pairwise
	global aligner algorithm and inherets from the PairwiseAligner
	class. This class can be used to align any two protein sequences.
	for NeedlemanWunsch

	Parameters
	----------
	substitution_matrix_file : str
		path/filename of substitution matrix for NeedlemanWunsch
    gap_open : float
		gap opening penalty for NeedlemanWunsch
    gap_extend : float
		gap extension penalty for NeedlemanWunsch

	Attributes
	----------
	seq1 : str
		sequence from fasta1 for NeedlemanWunsch
	seq2 : str
		sequence from fasta2 for NeedlemanWunsch
	header1 : str
		fasta header from fasta1 for NeedlemanWunsch
	header2 : str
		fasta header from fasta2 for NeedlemanWunsch
	alignment_ : tuple
		Tuple of format (seq1_alignment, seq2_alignment) for NeedlemanWunsch
	alignment_score_ : float
		Score of alignment from algorithm for NeedlemanWunsch

	"""
	def __init__(self, substitution_matrix_file, gap_open, gap_extend):
		super(NeedlemanWunsch, self).__init__(substitution_matrix_file, gap_open, gap_extend)

	def _fill_scoring_matrices(self, s_mat, x_mat, y_mat, tb_mat):
		"""This function uses dynamic programming to fill in all of the
		scoring and traceback matrices. This specific implementation
		is for global alignment using the NeedlemanWunsch alignment
		algorithm. The function takes in all initalized scoring matrices
		and fills them out. for SmithWaterman

		Parameters
		----------
		s_mat : array-like
			initialized core scoring matrix for NeedlemanWunsch
		x_mat : array-like
			initialized X gap scoring matrix, seq2 gapsm for NeedlemanWunsch
		y_mat : array-like
			initialized Y gap scoring matrix, seq1 gaps for NeedlemanWunsch
		tb_mat : array-like
			initalized traceback matrix

		Returns
		-------
		s_mat : array-like
			Filled core scoring matrix for NeedlemanWunsch
		x_mat : array-like
			Filled X gap scoring matrix, seq2 gaps for NeedlemanWunsch
		y_mat : array-like
			Filled Y gap scoring matrix, seq1 gaps for NeedlemanWunsch
		tb_mat : array-like
			Filled Traceback matrix, used for constructing alignment
			for NeedlemanWunsch

		"""
		# Filling in score Matrix
		seq1 = self.seq1
		seq2 = self.seq2
		for i in range(s_mat.shape[0])[1:]:
			for j in range(s_mat.shape[1])[1:]:
				# x gap matrix filling
				gap_dict = {1: 'new_gap', 0: 'extended_gap'}
				new_gap = (s_mat[i-1, j] - self.gap_open)
				extended_gap = (x_mat[i-1, j] - self.gap_extend)
				x_mat[i, j] = max([extended_gap, new_gap])
				gap_type_x = gap_dict[np.argmax([extended_gap, new_gap])]
				# y gap matrix filling
				new_gap = (s_mat[i, j-1] - self.gap_open)
				extended_gap = (y_mat[i, j-1] - self.gap_extend)
				y_mat[i, j] = max([extended_gap, new_gap])
				gap_type_y = gap_dict[np.argmax([extended_gap, new_gap])]
				# Taking maximum transition
				top = x_mat[i, j]
				side = y_mat[i, j]
				diag = s_mat[i-1, j-1] + self.sub_matrix[(seq1[i-1], seq2[j-1])]
				s_mat[i,j] = max((top, side, diag))
				max_ind = np.argmax((top, side, diag))
				# filling in traceback matrix
				if max_ind == 0 and gap_type_x == 'new_gap':
					tb_mat[i, j] = 'top_new'
				elif max_ind == 0 and gap_type_x == 'extended_gap':
					tb_mat[i, j] = 'top_extended'
				elif max_ind == 1 and gap_type_y == 'new_gap':
					tb_mat[i, j] = 'side_new'
				elif max_ind == 1 and gap_type_y == 'extended_gap':
					tb_mat[i, j] = 'side_extended'
				elif max_ind == 2:
					tb_mat[i, j] = 'diag'
		self.align_score_ = s_mat[-1, -1]
		return s_mat, x_mat, y_mat, tb_mat

	def _traceback(self, s_mat, tb_mat):
		"""This functions uses the core scoring matrix and the traceback
		matrix to find the alignment score and full alignment of seq1
		and seq2. This function is implemented for the global alignment
		NeedlemanWunsch algorithm. for NeedlemanWunsch

		Parameters
		----------
		s_mat : array-like
			Core scoring matrix for NeedlemanWunsch
		tb_mat : array-like
			Traceback matrix, used for constructing alignment for NeedlemanWunsch

		Returns
		-------
		alignment : tuple
			Tuple of format (seq1_alignment, seq2_alignment) for NeedlemanWunsch
		alignment_score : float
			Float value of maximum alignment score in the core score matrix
			for NeedlemanWunsch

		"""
		# Traceback
		seq1 = self.seq1
		seq2 = self.seq2
		alignment_seq1 = ''
		alignment_seq2 = ''
		i = len(seq1)
		j = len(seq2)
		while i > 0 and j > 0:
			tb_val = tb_mat[i, j]
			if tb_val == 'top_extended':
				alignment_seq2 = '-' + alignment_seq2
				alignment_seq1 = seq1[i-1] + alignment_seq1
				i -= 1
			elif tb_val == 'top_new':
				alignment_seq2 = '-' + alignment_seq2
				alignment_seq1 = seq1[i-1] + alignment_seq1
				i -= 1
			elif tb_val == 'side_extended':
				alignment_seq2 = seq2[j-1] + alignment_seq2
				alignment_seq1 = '-' + alignment_seq1
				j -= 1
			elif tb_val == 'side_new':
				alignment_seq2 = seq2[j-1] + alignment_seq2
				alignment_seq1 = '-' + alignment_seq1
				j -= 1
			elif tb_val == 'diag':
				alignment_seq2 = seq2[j-1] + alignment_seq2
				alignment_seq1 = seq1[i-1] + alignment_seq1
				i -= 1
				j -= 1
		self.alignment_ = (alignment_seq1, alignment_seq2)
		alignment = self.alignment_
		alignment_score = self.align_score_
		return alignment, alignment_score

	def _algorithm_matrix_initiation(self, s_mat):
		"""This function intialized the core score matrix for the
		NeedlemanWunsch Alignment by setting the first row and
		column to be set to an approximation of negative
		infinitiy. for NeedlemanWunsch

		Parameters
		----------
		s_mat : array-like
			Initialized core scoring matrix for NeedlemanWunsch

		Returns
		-------
		s_mat : array-like
			Core scoring matrix initalized for NeedlemanWunsch

		"""
		# Filling in first column and row
		seq1 = self.seq1
		seq2 = self.seq2
		s_mat[1:, 0] = np.arange(len(seq1)+1)[:-1]*-self.gap_open*100000000
		s_mat[0, 1:] = np.arange(len(seq2)+1)[:-1]*-self.gap_open*100000000
		return s_mat
