##### UPDATES THAT NEED TO HAPPEN #####
# 1) Modify traceback to allowing branching (ask others)
# 2) Modify align returned attributes to allow branching
# 3) Write class doc strings
# 4) Write up answers to questions in Readme.md
# 5

class PairwiseAligner:
	'''
	This is the base class for all pairwise aligment algorithms
	implemented in this API. All pariwise alignment classes
	inheret from this base class.
	'''
	def __init__(self, substitution_matrix_file):
		# generating substitution matrix
		self.sub_matrix = _initiate_sub_matrix(substitution_matrix_file)
		self.seq1 = None
		self.seq2 = None
		self.score_matrix = None

	def _initate_sub_matrix(self, matrix_file):
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
		with open(file, 'r') as f:
			dict_sm = {} # Dictionary for storing scores
			residue_list = [] # For storing residue list
			start = False # trigger for reading in score values
			res_2 = 0 # used for generating substitution matrix
			# reading file line by line
			for count, line in enumerate(f):
	            # Reading in Residue List
				if line.strip().startswith('A ') and start == False:
					for val in line.strip():
						if val != ' ':
							residue_list.append(val)
							# if residue line encountered then read in matrix
							start = True
							# Generating substitution scoring dictionary
						elif start == True and res_2 < len(residue_list):
							line = [k for k in line.strip().split(' ') if k != '']
							line.strip().replace('  ', ' ').split(' ')
							# reading in line by line to create substitution dictionary
							for res_1 in range(len(line))[res_2:]:
								dict_sm[(residue_list[res_1], residue_list[res_2])] = line[res_1]
								dict_sm[(residue_list[res_2], residue_list[res_1])] = line[res_1]
								res_2 += 1
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
			for line in file:
				if is_header and line.strip().startswith('>'):
					header = line.strip() # reading in fasta header
					is_header = False
				else:
					seq += line.strip() # generating full seq
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

	def _initiate_score_matrix(self):
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
		seq1_len = len(self.seq1)
		seq2_len = len(self.seq2)
		score_matrix = np.zeros((seq1_len, seq2_len))
		return score_matrix

	def align(self):
		'''
		This is a method that must be overwritten when inhereting from
		this class.
		'''
		pass


class SmithWaterman(PairwiseAligner):
	def __init__(self, substitution_matrix_file, gap_open, gap_extend):
		super(SmithWaterman, self).__init__(substitution_matrix_file)
		self.gap_open = gap_open
		self.gap_extend = gap_extend

	def align(self, fasta1, fasta2):
		'''
		This funciton implements the dynamic programming
		and actaul local alignment process of SmithWaterman
		alignmnet. The alignment score is saved to the field
		self.align_score_ and the alignments are saved as
		a list of composition [seq1_alignment, seq2_alignment]
		to the field self.aligment_.

		Args
		----
		fasta1 (str): Filename and path to first fasta file for alignment.
		fasta2 (str): Filename and path to second fasta file for alignment.

		Returns
		-------
		None
		'''
		# Initiating Sequences and Scoring Matirx
		self._initiate_seqs(fasta_1, fasta_2)
		self.score_matrix = self._initiate_score_matrix()
		s_matrix = self.score_matrix
		t_matrix = np.zeros(score_matrix.shape) # seq1 gap matrix
		u_matrix = np.zeros(score_matrix.shape) # seq2 gap matrix
		# filling out score matrix (forward pass)
		for i in range(s_matrix.shape[0])[1:]:
			for j in range(s_matrix.shape[1])[1:]:
				# t gap matrix filling
				t_matrix[i, j] = max([(t_matrix[i-1, j] - self.gap_extend), (s_matrix[i-1, j] - (self.gap_extend + self.gap_open)), 0])
				# u gap matrix filling
				u_matrix[i, j] = max([(u_matrix[i, j-1] - self.gap_extend), (s_matrix[i, j-1] - (self.gap_extend + self.gap_open)), 0])
				top = t[i, j]
				side = u[i, j]
				diag = s[i-1, j-1] + self.sub_matrix[(seq1[i], seq2[j])]
				s_matrix[i,j] = max((top, side, diag, 0))
		# Traceback (backwards pass)
		max_ind = np.unravel_index(np.argmax(s_matrix, axis=None), s_matrix.shape)
		self.align_score_ = np.amax(s_matrix)
		i, j = max_ind
		alignment_seq1 = seq1[i]
		alignment_seq2 = seq2[j]
		while i > 0 and j > 0:
			top = s_matrix[i-1, j]
			side = s_matrix[i, j-1]
			diag = s_matrix[i-1, j-1]
			best_ind = np.argmax(np.array([top, side, diag]))
			best_max = np.amax(np.array([top, side, diag]))
			if best_max == 0:
				break
			if best_ind == 0:
				i -= 1
				alignment_seq1 += seq1[i]
				alignment_seq2 += '-'
			elif best_ind == 1:
				j -= 1
				alignment_seq1 += '-'
				alignment_seq2 += seq2[j]
			elif best_ind == 2:
				i -= 1
				j -= 1
				alignment_seq1 += seq1[i]
				alignment_seq2 += seq2[j]
		self.alignment_ = [alignment_seq1, alignment_seq2]


class NeedlemanWunsch(PairwiseAligner):
	def __init__(self, substitution_matrix_file, gap_penalty):
		super(NeedlemanWunsch, self).__init__(substitution_matrix_file)
		self.gap_penalty = gap_penalty

	def align(self, fasta1, fasta2):
		'''
		This funciton implements the dynamic programming
		and actaul global alignment process of Needleman Wunch
		alignmnet. The alignment score is saved to the field
		self.align_score_ and the alignments are saved as
		a list of composition [seq1_alignment, seq2_alignment]
		to the field self.aligment_.

		Args
		----
		fasta1 (str): Filename and path to first fasta file for alignment.
		fasta2 (str): Filename and path to second fasta file for alignment.

		Returns
		-------
		None
		'''
		# Initiating Sequences and Scoring Matirx
		self._initiate_seqs(fasta_1, fasta_2)
		self.score_matrix = self._initiate_score_matrix()
		s_matrix = self.score_matrix
		# Filling in score Matrix
		for i in range(s_matrix.shape[0]):
			for j in range(s_matrix.shape[1]):
				if i == j and i == 0:
					pass
				else:
					top = s_matrix[i-1, j] - self.gap_penalty
					side = s_matrix[i, j-1] - self.gap_penalty
					diag = s_matrix[i-1, j-1] + self.sub_matrix[(seq1[i], seq2[j])]
					s_matrix[i, j] = max([top, side, diag])
		# Traceback
		alignment_seq1 = seq1[-1]
		alignment_seq2 = seq2[-1]
		###### NEED TO MODIFY THIS ######
		# max_ind = np.unravel_index(np.argmax(s_matrix, axis=None), s_matrix.shape)
		# self.align_score_ = np.amax(s_matrix)
		# i, j = max_ind
		# alignment_seq1 = seq1[i]
		# alignment_seq2 = seq2[j]
		while i > 0 and j > 0:
			top = s[i-1, j]
			side = s[i, j-1]
			diag = s[i-1, j-1]
			best_ind = np.argmax(np.array([top, side, diag]))
			best_max = np.amax(np.array([top, side, diag]))
			if best_ind == 0:
				i -= 1
				alignment_seq1 += seq1[i]
				alignment_seq2 += '-'
			elif best_ind == 1:
				j -= 1
				alignment_seq1 += '-'
				alignment_seq2 += seq2[j]
			elif best_ind == 2:
				i -= 1
				j -= 1
				alignment_seq1 += seq1[i]
				alignment_seq2 += seq2[j]
		self.alignment_ = [alignment_seq1, alignment_seq2]
