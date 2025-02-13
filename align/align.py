# Importing Dependencies
import numpy as np
from typing import Tuple

# Defining class for Needleman-Wunsch Algorithm for Global pairwise alignment
class NeedlemanWunsch:
    """ Class for NeedlemanWunsch Alignment

    Parameters:
        sub_matrix_file: str
            Path/filename of substitution matrix
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty

    Attributes:
        seqA_align: str
            seqA alignment
        seqB_align: str
            seqB alignment
        alignment_score: float
            Score of alignment from algorithm
        gap_open: float
            Gap opening penalty
        gap_extend: float
            Gap extension penalty
    """
    def __init__(self, sub_matrix_file: str, gap_open: float, gap_extend: float):
        # Init alignment and gap matrices
        self._align_matrix = None
        self._gapA_matrix = None
        self._gapB_matrix = None

        # Init matrices for backtrace procedure
        self._back = None
        self._back_A = None
        self._back_B = None

        # Init alignment_score
        self.alignment_score = 0

        # Init empty alignment attributes
        self.seqA_align = ""
        self.seqB_align = ""

        # Init empty sequences
        self._seqA = ""
        self._seqB = ""

        # Setting gap open and gap extension penalties
        self.gap_open = gap_open
        assert gap_open < 0, "Gap opening penalty must be negative."
        self.gap_extend = gap_extend
        assert gap_extend < 0, "Gap extension penalty must be negative."

        # Generating substitution matrix
        self.sub_dict = self._read_sub_matrix(sub_matrix_file) # substitution dictionary

    def _read_sub_matrix(self, sub_matrix_file):
        """
        DO NOT MODIFY THIS METHOD! IT IS ALREADY COMPLETE!

        This function reads in a scoring matrix from any matrix like file.
        Where there is a line of the residues followed by substitution matrix.
        This file also saves the alphabet list attribute.

        Parameters:
            sub_matrix_file: str
                Name (and associated path if not in current working directory)
                of the matrix file that contains the scoring matrix.

        Returns:
            dict_sub: dict
                Substitution matrix dictionary with tuple of the two residues as
                the key and score as value e.g. {('A', 'A'): 4} or {('A', 'D'): -8}
        """
        with open(sub_matrix_file, 'r') as f:
            dict_sub = {}  # Dictionary for storing scores from sub matrix
            residue_list = []  # For storing residue list
            start = False  # trigger for reading in score values
            res_2 = 0  # used for generating substitution matrix
            # reading file line by line
            for line_num, line in enumerate(f):
                # Reading in residue list
                if '#' not in line.strip() and start is False:
                    residue_list = [k for k in line.strip().upper().split(' ') if k != '']
                    start = True
                # Generating substitution scoring dictionary
                elif start is True and res_2 < len(residue_list):
                    line = [k for k in line.strip().split(' ') if k != '']
                    # reading in line by line to create substitution dictionary
                    assert len(residue_list) == len(line), "Score line should be same length as residue list"
                    for res_1 in range(len(line)):
                        dict_sub[(residue_list[res_1], residue_list[res_2])] = float(line[res_1])
                    res_2 += 1
                elif start is True and res_2 == len(residue_list):
                    break
        return dict_sub
    
    def _initialize_matrix(self, m, n):
        '''
        Intializing all relevant matrices
        
        '''
        self._align_matrix = np.zeros((n+1, m+1)) # initialize the alignment matrix with score of 0
        self._gapA_matrix = np.zeros((n+1, m+1))
        self._gapB_matrix = np.zeros((n+1, m+1))

        # initializing the first column/row based on the gap open and extend penalties
        for i in range(n+1):
            # Initializing the first column of the gapA matrix with -inf and the gapB matrix with corresponding scores
            self._gapA_matrix[i][0] = self.gap_open + i*self.gap_extend
            self._gapB_matrix[i][0] = -np.inf

            # Initializing the alignment matrix with corresponding score
            if i <= 1:
                self._align_matrix[i][0] = i*self.gap_open
            else:
                self._align_matrix[i][0] = self.gap_open + (i-1)*self.gap_extend
        

        for j in range(m+1):
            # Initializing gap matrices (top and left edges, making sure not to over-write the top corner for gapA matrix (set above))
            if j != 0:
                self._gapA_matrix[0][j] = -np.inf
            self._gapB_matrix[0][j] = self.gap_open + j*self.gap_extend

            if j <= 1:
                self._align_matrix[0][j] = j*self.gap_open
            else:
                self._align_matrix[0][j] = self.gap_open + (j-1)*self.gap_extend       


        # Initialize backtrace matrix
        self._back = np.full((n+1, m+1), 'NaN', dtype=str)
        self._back_A = np.full((n+1, m+1), 'NaN', dtype=str)
        self._back_B = np.full((n+1, m+1), 'NaN', dtype=str)


    def align(self, seqA: str, seqB: str) -> Tuple[float, str, str]:
        """
        TODO
        
        This function performs global sequence alignment of two strings
        using the Needleman-Wunsch Algorithm
        
        Parameters:
        	seqA: str
         		the first string to be aligned
         	seqB: str
         		the second string to be aligned with seqA
         
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
                 

        This function is implemented with guidance of the YouTube video provided in the README
        """
        # Resetting alignment in case method is called more than once
        self.seqA_align = ""
        self.seqB_align = ""

        # Resetting alignment score in case method is called more than once
        self.alignment_score = 0

        # Initializing sequences for use in backtrace method
        self._seqA = seqA.upper()
        self._seqB = seqB.upper()
        
        # TODO: Initialize matrix private attributes for use in alignment
        # create matrices for alignment scores, gaps, and backtracing
        
        # matrix will by n + 1 x m + 1 where n is the length of seq A and m is the lenght of seq b
        # n+1 x m+1 since we need the first column/row to be a "gap"
        # seqA will be across (first index) and seqB will be down (second index)

        m, n = len(self._seqA), len(self._seqB)
        self._seqA_len = m
        self._seqB_len = n

        self._initialize_matrix(m=m, n=n)

        # filling rest of alignment score matrix
        # M(i,j) = max(
        #               M(i-1, j) + penalty, # Seq A aligns with gap of seq B
        #               M(i, j-1) + penalty, # seq b aligns with gap of seq A
        #               M(i-1, j-1) + match/mismatch # seq a + b aligns with each other (could still be mismatch)
        #               )

        # Iterating over the whole table to create alignment score tables + backtrace tables
        for A_ix in range(1, m+1):
            A_aa = self._seqA[A_ix-1]
            for B_ix in range(1, n+1):
                B_aa = self._seqB[B_ix-1]
                match_score = self.sub_dict[(A_aa, B_aa)]

                # M scores
                M_diag_score = self._align_matrix[B_ix-1][A_ix-1] + match_score
                M_a_gap_score = self._gapA_matrix[B_ix-1][A_ix-1] + match_score
                M_b_gap_score = self._gapB_matrix[B_ix-1][A_ix-1] + match_score

                M_scores = [M_diag_score, M_a_gap_score, M_b_gap_score]
                M_max_ix = np.argmax(M_scores)

                self._align_matrix[B_ix][A_ix] = M_scores[M_max_ix]
                
                # Gap A matrix scores
                A_gap_open_score = self._align_matrix[B_ix-1][A_ix] + self.gap_open + self.gap_extend
                A_gap_extend_score = self._gapA_matrix[B_ix-1][A_ix] + self.gap_extend

                A_scores = [A_gap_open_score, A_gap_extend_score]
                A_max_ix = np.argmax(A_scores)

                self._gapA_matrix[B_ix][A_ix] = A_scores[A_max_ix]

                # Gap B matrix scores
                B_gap_open_score = self._align_matrix[B_ix][A_ix-1] + self.gap_open + self.gap_extend
                B_gap_extend_score = self._gapB_matrix[B_ix][A_ix-1] + self.gap_extend

                B_scores = [B_gap_open_score, B_gap_extend_score]
                B_max_ix = np.argmax(B_scores)

                self._gapB_matrix[B_ix][A_ix] = B_scores[B_max_ix]

                '''
                To make backtrace table:

                for M_max_ix: if 0-> came from align_M, if 1-> came from gapA, 2 -> came from gapB (always comes from diagonal)
                for gapA: if 0 -> came from align_M, if 1 -> came from gapA (always from above)
                for gapB: if 0 -> came from align_M, if 1-> came from gapB (always from left)
                '''
                if A_max_ix == 0:
                    a_backtrace = "M"
                elif A_max_ix == 1:
                    a_backtrace = "A"

                if B_max_ix == 0:
                    b_backtrace = "M"
                elif B_max_ix == 1:
                    b_backtrace = "B"

                if M_max_ix == 0:
                    m_backtrace = "M"
                elif M_max_ix == 1:
                    m_backtrace = "A"
                elif M_max_ix == 2:
                    m_backtrace = "B"

                self._back[B_ix][A_ix] = m_backtrace
                self._back_A[B_ix][A_ix] = a_backtrace
                self._back_B[B_ix][A_ix] = b_backtrace
                
        self._align_score = np.max([self._align_matrix[n][m], self._gapA_matrix[n][m], self._gapB_matrix[n][m]])
        		    
        return self._backtrace()

    def _backtrace(self) -> Tuple[float, str, str]:
        """
        TODO
        
        This function traces back through the back matrix created with the
        align function in order to return the final alignment score and strings.
        
        Parameters:
        	None
        
        Returns:
         	(alignment score, seqA alignment, seqB alignment) : Tuple[float, str, str]
         		the score and corresponding strings for the alignment of seqA and seqB
        """
        m = self._seqA_len
        n = self._seqB_len

        m_ix = m
        n_ix = n
        index_pair_list = []
        prev_matrix = ''

        while m_ix > 0 or n_ix > 0:
            m_score = self._align_matrix[n_ix][m_ix]
            a_score = self._gapA_matrix[n_ix][m_ix]
            b_score = self._gapB_matrix[n_ix][m_ix]

            score_list = [m_score, a_score, b_score]
            max_ix = np.argmax(score_list)

            index_pair_list.append((m_ix, n_ix))

            if max_ix == 0:
                # comes from align matrix
                prev_matrix = self._back[n_ix][m_ix]

                m_ix -= 1
                n_ix -= 1

                seqA_aa = self._seqA[m_ix]
                seqB_aa = self._seqB[n_ix]

            elif max_ix == 1:
                prev_matrix = self._back_A[n_ix][m_ix]
                n_ix -= 1

                seqA_aa = "-"
                seqB_aa = self._seqB[n_ix]


            elif max_ix == 2:
                prev_matrix = self._back_B[n_ix][m_ix]
                m_ix -= 1

                seqA_aa = self._seqA[m_ix]
                seqB_aa = '-'

            if m_ix < 0:
                seqA_aa = "-"
            if n_ix < 0:
                seqB_aa = "-"

            self.seqA_align += seqA_aa
            self.seqB_align += seqB_aa

        self.seqA_align = self.seqA_align[::-1]
        self.seqB_align = self.seqB_align[::-1]
        self.alignment_score = self._align_score

        return (self.alignment_score, self.seqA_align, self.seqB_align)


def read_fasta(fasta_file: str) -> Tuple[str, str]:
    """
    DO NOT MODIFY THIS FUNCTION! IT IS ALREADY COMPLETE!

    This function reads in a FASTA file and returns the associated
    string of characters (residues or nucleotides) and the header.
    This function assumes a single protein or nucleotide sequence
    per fasta file and will only read in the first sequence in the
    file if multiple are provided.

    Parameters:
        fasta_file: str
            name (and associated path if not in current working directory)
            of the Fasta file.

    Returns:
        seq: str
            String of characters from FASTA file
        header: str
            Fasta header
    """
    assert fasta_file.endswith(".fa"), "Fasta file must be a fasta file with the suffix .fa"
    with open(fasta_file) as f:
        seq = ""  # initializing sequence
        first_header = True
        for line in f:
            is_header = line.strip().startswith(">")
            # Reading in the first header
            if is_header and first_header:
                header = line.strip()  # reading in fasta header
                first_header = False
            # Reading in the sequence line by line
            elif not is_header:
                seq += line.strip().upper()  # generating full sequence
            # Breaking if more than one header is provided in the fasta file
            elif is_header and not first_header:
                break
    return seq, header