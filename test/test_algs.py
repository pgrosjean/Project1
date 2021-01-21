import pytest
from align import algs

@pytest.fixture
def some_relevant_data():
	return np.ones(10)

def test_fasta_io_sw():
	sw_aligner = algs.SmithWaterman('../scoring_matrices/BLOSUM50.mat', [4, 1])
	seq_t1, header_t1 = sw_aligner._read_fasta('test1.fa')
	seq_t2, header_t2 = sw_aligner._read_fasta('test2.fa')
	assert seq_t1 == 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'
	assert seq_t2 == 'GAWDDAMDAWG'
	assert header_t1 == '>test_1_fasta'
	assert header_t2 == '>test_2_fasta'

def test_fasta_io_nw():
	nw_aligner = algs.NeedlemanWunsch('../scoring_matrices/BLOSUM50.mat', 1)
	seq_t1, header_t1 = nw_aligner._read_fasta('test1.fa')
	seq_t2, header_t2 = nw_aligner._read_fasta('test2.fa')
	assert seq_t1 == 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'
	assert seq_t2 == 'GAWDDAMDAWG'
	assert header_t1 == '>test_1_fasta'
	assert header_t2 == '>test_2_fasta'

def test_seq_initiation_sw():
	sw_aligner = algs.SmithWaterman('../scoring_matrices/BLOSUM50.mat', [4, 1])
	sw_aligner._initiate_seqs('test1.fa', 'test2.fa')
	assert sw_aligner.seq1 == 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'
	assert sw_aligner.seq2 == 'GAWDDAMDAWG'
	assert sw_aligner.header1 == '>test_1_fasta'
	assert sw_aligner.header2 == '>test_2_fasta'

def test_seq_initiation_nw():
	nw_aligner = algs.NeedlemanWunsch('../scoring_matrices/BLOSUM50.mat', 1)
	nw_aligner._initiate_seqs('test1.fa', 'test2.fa')
	assert nw_aligner.seq1 == 'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM'
	assert nw_aligner.seq2 == 'GAWDDAMDAWG'
	assert nw_aligner.header1 == '>test_1_fasta'
	assert nw_aligner.header2 == '>test_2_fasta'

def test_score_matrix_gen_sw():
	sw_aligner = algs.SmithWaterman('../scoring_matrices/BLOSUM50.mat', [4, 1])
	sw_aligner._initiate_seqs('test1.fa', 'test2.fa')
	sw_aligner.score_matrix = self._initiate_score_matrix()
	seq1 = nw_aligner.seq1
	seq2 = nw_alginer.seq2
	assert sw_aligner.score_matrix.shape == (len(seq1), len(seq2))

def test_score_matrix_gen_nw():
	nw_aligner = algs.NeedlemanWunsch('../scoring_matrices/BLOSUM50.mat', 1)
	nw_aligner._initiate_seqs('test1.fa', 'test2.fa')
	nw_aligner.score_matrix = self._initiate_score_matrix()
	seq1 = nw_aligner.seq1
	seq2 = nw_alginer.seq2
	assert nw_aligner.score_matrix.shape == (len(seq1), len(seq2))

def test_scoring_matrix_io_sw():
	assert True

def test_scoring_matrix_io_nw():
	assert True

def test_identical_sw():
	sw_aligner = algs.SmithWaterman('../scoring_matrices/BLOSUM62.mat', [10, 0.5])
	sw_aligner.align('test3.fa', 'test3.fa')
	assert sw_aligner.align_score_ == 733.0
	assert sw_aligner.alignment_ == ['SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM',
	'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM']
	assert True

def test_identical_nw():
	nw_alginer = algs.NeedlemanWunsch('../scoring_matrices/BLOSUM62.mat', 1)
	nw_alginer.align('test3.fa', 'test3.fa')
	assert nw_aligner.align_score_ == 734.0
	assert nw_aligner.alignment_ == ['SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM',
	'SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM']

def test_alignment_score_sw():
	sw_aligner = algs.SmithWaterman('../scoring_matrices/BLOSUM62.mat', [10, 0.5])
	sw_aligner.align('test3.fa', 'test4.fa')
	assert sw_aligner.align_score_ == 648.0
	assert sw_aligner.alignment_ == ['SLEAAQKSNVTSSWAKASAAWGTAGPEFFMALFDAHDDVFAKFSGLFSGAAKGTVKNTPEMAAQAQSFKGLVSNWVDNLDNAGALEGQCKTFAANHKARGISAGQLEAAFKVLSGFMKSYGGDEGAWTAVAGALMGEIEPDM',
	'MVLSGEDKSNIKAAWGKIGGHGAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHGKKVADALASAAGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASHHPADFTPAVHASLDKFLASVSTVLTSKYR']

def test_alignment_score_nw():
	nw_aligner = algs.NeedlemanWunsch('../scoring_matrices/BLOSUM62.mat', 1)
	nw_aligner.align('test3.fa', 'test4.fa')
	assert sw_aligner.align_score_ == 656.0
	assert sw_aligner.alignment_ == ['MVLSPA-DKTNVKAAWGKVGAHAG-EYGAEALERMFLSFPTTKTYFPHFDLSHGSAQVKGHGKKVADALTNAVA-HVDDMPNALSALSDLHAHKLRVDPVNFKLLSHCLLVTLAAHL-PAEFTPAVHASLDKFLASVSTVLTSKYR',
	'MVLS-GEDKSNIKAAWGKIGGH-GAEYGAEALERMFASFPTTKTYFPHFDVSHGSAQVKGHGKKVADALASA-AGHLDDLPGALSALSDLHAHKLRVDPVNFKLLSHCLLVTLASH-HPADFTPAVHASLDKFLASVSTVLTSKYR']
