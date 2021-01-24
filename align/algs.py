##### UPDATES THAT NEED TO HAPPEN #####
# Implement backtrack table for all of these

import numpy as np

class PairwiseAligner:
	'''
	This is the base class for all pairwise alignment algorithms
	implemented in this API. All pariwise alignment classes
	inheret from this base class.
	'''
	def __init__(self, substitution_matrix_file, gap_open, gap_extend):
		# generating substitution matrix
		self.sub_matrix = self._initiate_sub_matrix(substitution_matrix_file)
		self.gap_open = gap_open
		self.gap_extend = gap_extend
		self.seq1 = None
		self.seq2 = None
		self.score_matrix = None

	def _initiate_sub_matrix(self, matrix_file):
		'''
		This function reads in a scoring matrix from any matrix like file.
		Where there is a line of the residues followed by substitution
		matrix.

		Args
		----
		matrix_file (str): name (and associated path if not in current
			working directory) of the matrix file that contains the
			scoring matrix.

		Returns
		-------
		dict_sm (dict): substitution matrix dictionary with tuple
			of the two residues as the key and score as value
			e.g. {('A', 'A') : 4} or {('A', 'D') : -8}
		'''
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
		'''
		This function reads in a FASTA file and returns the associated
		string of characters (residues or nucleotides) and the header.

		Args
		----
		fasta_file (str): name (and associated path if not in current
			working directory) of the Fasta file.

		Returns
		-------
		seq (str): String of characters from FASTA file
		header (str): Fasta header
		'''
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
		'''
		This function initiates the sequences for alignment by
		calling the _read_fasta method on the two fasta files.
		This saves the sequences and headers into the fields
		self.seq1, self.seq2, self.header1, self.header2.

		Args
		----
		fasta_1 (str): path and filename of first fasta for alignment.
		fasta_2 (str): path and filename of second fasta for alignment.

		Returns
		-------
		None
		'''
		seq1, header1 = self._read_fasta(fasta_1)
		seq2, header2 = self._read_fasta(fasta_2)
		self.seq1 = seq1
		self.seq2 = seq2
		self.header1 = header1
		self.header2 = header2

	def initiate_matrix(self):
		'''
		This function initializes the score matrix of dimensions
		len(seq1)xlen(seq2). This method must be run after the
		_initate_seqs method.

		Args
		----
		None

		Returns
		-------
		score_matrix (numpy array): array of size len(seq1)xlen(seq2)
		'''
		seq1_len = len(self.seq1) + 1
		seq2_len = len(self.seq2) + 1
		score_matrix = np.zeros((seq1_len, seq2_len))
		return score_matrix

	def align(self):
		'''
		This is a method that must be overwritten when inhereting from
		this class.
		'''
		pass



################################################################################




class SmithWaterman(PairwiseAligner):
	def __init__(self, substitution_matrix_file, gap_open, gap_extend):
		super(SmithWaterman, self).__init__(substitution_matrix_file, gap_open, gap_extend)

	def align(self, fasta1, fasta2):
		'''
		This funciton implements the dynamic programming
		and actaul local alignment process of SmithWaterman
		alignmnet. The alignment score is saved to the field
		self.align_score_ and the alignments are saved as
		a list of composition [seq1_alignment, seq2_alignment]
		to the field self.alignment_.

		Args
		----
		fasta1 (str): Filename and path to first fasta file for alignment.
		fasta2 (str): Filename and path to second fasta file for alignment.

		Returns
		-------
		None
		'''
		# Initializing Sequences
		self._initiate_seqs(fasta1, fasta2)
		seq1 = self.seq1
		seq2 = self.seq2
		# Initializing Core Score Matrix
		self.score_matrix = self.initiate_matrix()
		s_mat = self.score_matrix
		# Initializing Additional Scoring Matrices for Indels (Gaps)
		x_mat = self.initiate_matrix() # seq2 gap matrix
		x_mat[:, 0] = np.arange(len(seq1)+1)*-self.gap_extend - self.gap_open
		x_mat[0, :] = np.arange(len(seq2)+1)*-self.gap_open*100000000
		y_mat = self.initiate_matrix() # seq1 gap matrix
		y_mat[:, 0] = np.arange(len(seq1)+1)*-self.gap_open*100000000
		y_mat[0, :] = np.arange(len(seq2)+1)*-self.gap_extend - self.gap_open
		# Initialize Traceback matrix
		tb_mat = self.initiate_matrix().astype('str')
		# Filling in first column and row
		s_mat[1:, 0] = 0
		s_mat[0, 1:] = 0

		# Filling in score Matrix
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
		# Traceback
		alignment_seq1 = ''
		alignment_seq2 = ''
		self.align_score_ = s_mat[-1, -1]
		max_ind = np.unravel_index(np.argmax(s_mat, axis=None), s_mat.shape)
		self.align_score_ = np.amax(s_mat)
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
		self.alignment_ = [alignment_seq1, alignment_seq2]


class NeedlemanWunsch(PairwiseAligner):
	def __init__(self, substitution_matrix_file, gap_open, gap_extend):
		super(NeedlemanWunsch, self).__init__(substitution_matrix_file, gap_open, gap_extend)

	def align(self, fasta1, fasta2):
		'''
		This funciton implements the dynamic programming
		and actaul global alignment process of Needleman Wunch
		alignmnet. The alignment score is saved to the field
		self.align_score_ and the alignments are saved as
		a list of composition [seq1_alignment, seq2_alignment]
		to the field self.alignment_.

		Args
		----
		fasta1 (str): Filename and path to first fasta file for alignment.
		fasta2 (str): Filename and path to second fasta file for alignment.

		Returns
		-------
		None
		'''
		# Initializing Sequences
		self._initiate_seqs(fasta1, fasta2)
		seq1 = self.seq1
		seq2 = self.seq2
		# Initializing Core Score Matrix
		self.score_matrix = self.initiate_matrix()
		s_mat = self.score_matrix
		# Initializing Additional Scoring Matrices for Indels (Gaps)
		x_mat = self.initiate_matrix() # seq2 gap matrix
		x_mat[:, 0] = np.arange(len(seq1)+1)*-self.gap_extend - self.gap_open
		x_mat[0, :] = np.arange(len(seq2)+1)*-self.gap_open*100000000
		y_mat = self.initiate_matrix() # seq1 gap matrix
		y_mat[:, 0] = np.arange(len(seq1)+1)*-self.gap_open*100000000
		y_mat[0, :] = np.arange(len(seq2)+1)*-self.gap_extend - self.gap_open
		# Initialize Traceback matrix
		tb_mat = self.initiate_matrix().astype('str')
		# Filling in first column and row
		s_mat[1:, 0] = np.arange(len(seq1)+1)[:-1]*-self.gap_open*100000000
		s_mat[0, 1:] = np.arange(len(seq2)+1)[:-1]*-self.gap_open*100000000
		# Initializing traceback matrix
		tb_mat = self.initiate_matrix().astype('str')
		# Filling in score Matrix
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
		# Traceback
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
		self.alignment_ = [alignment_seq1, alignment_seq2]
		self.align_score_ = s_mat[-1, -1]
