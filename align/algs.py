
class PairwiseAligner:
	'''
	This is the base class for all pairwise aligment algorithms
	implemented in this API. All pariwise alignment classes
	inheret from this base class.
	'''
	def __init__(self, substitution_matrix_file, fasta_1, fasta_2):
		# generating substitution matrix
		self.sub_matrix = _initiate_sub_matrix(substitution_matrix_file)
		self._initiate_seqs(fasta_1, fasta_2) # initiates seq and header fields
		self.score_matrix = self._initiate_score_matrix()

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
                    line = line.strip().replace('  ', ' ').split(' ')
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
				if is_header:
					header = line.strip() # reading in fasta header
					is_header = False
				else:
					seq += line.strip() # generating full seq
		return seq, header

	def _initiate_seqs(self, fasta_1, fasta_2):
		'''
		This function initiates the sequences for alignment
		'''
		seq1, header1 = self._read_fasta(fasta_1)
		seq2, header2 = self._read_fasta(fasta_2)
		self.seq1 = seq1
		self.seq2 = seq2
		self.header1 = header1
		self.header2 = header2
		return None

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
		This is a method that must be overwritten
		when inhereting from this class.
		'''
		pass


class SmithWaterman(PairwiseAligner):
	def __init__(self, substitution_matrix_file, fasta_1, fasta_2, scoring_scheme):
		super(SmithWaterman).__init__(substitution_matrix_file, fasta_1, fasta_2)
		self.gap_open = scoring_scheme[0]
		self.gap_extend = scoring_scheme[1]

	def align(self):
		# creating local score_matrix
		s_matrix = self.score_matrix
		t_matrix = np.zeros(score_matrix.shape) # seq1 gap matrix
		u_matrix = np.zeros(score_matrix.shape) # seq2 gap matrix
		# filling out score matrix (forward pass)
		for i in range(s_matrix.shape[0])[1:]:
			for j in range(s_matrix.shape[1])[1:]:
				# t gap matrix filling
				t_matrix[i, j] = max([(t_matrix[i-1, j] + self.gap_extend), (s_matrix[i-1, j] - (self.gap_extend + self.gap_open))])
				# u gap matrix filling
				u_matrix[i, j] = max([(u_matrix[i, j-1] + self.gap_extend), (s_matrix[i, j-1] - (self.gap_extend + self.gap_open))])
				top = t[i, j]
				side = u[i, j]
				diag = s[i-1, j-1] + self.sub_matix[(seq1[i], seq2[j])]
				s_matrix[i,j] = max((top, side, diag))
		# traceback (backwards pass)
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
			best_max = np.max(np.array([top, side, diag]))
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
		self.alignment = [alignment_seq1, alignment_seq2]
		return None


class NeedlemanWunsch(PairwiseAligner):
	pass
