#from helper_functions import global_alignment
#pip install biopython matplotlib jupyterlab

def global_alignment(seq1, seq2, scoring_function):
    """Global sequence alignment using the Needleman–Wunsch algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> global_alignment("abracadabra", "dabarakadara", lambda x, y: [-1, 1][x == y])
    ('-ab-racadabra', 'dabarakada-ra', 5.0)

    Other alignments are not possible.

    """
   
    # Load BLOSUM62 substitution matrix
    blosum62 = MatrixInfo.blosum62

    n, m = len(seq1), len(seq2)

    # INITIALISATION
    score = [[0] * (m + 1) for _ in range(n + 1)]
    back = [[None] * (m + 1) for _ in range(n + 1)]
    gap_penalty = scoring_function('-', '-')

    for i in range(1, n + 1):
        score[i][0] = score[i - 1][0] + gap_penalty
        back[i][0] = "up"

    for j in range(1, m + 1):
        score[0][j] = score[0][j - 1] + gap_penalty
        back[0][j] = "left"

    # RECURRENCE
    for i in range(1, n + 1):
        for j in range(1, m + 1):

            a1 = seq1[i - 1]
            a2 = seq2[j - 1]

            # BLOSUM62 lookup
            match_score = blosum62.get((a1, a2))
            if match_score is None:
                match_score = blosum62.get((a2, a1))
            if match_score is None:
                match_score = -1  # fallback for unknown characters

            diag = score[i - 1][j - 1] + match_score
            up   = score[i - 1][j] + gap_penalty
            left = score[i][j - 1] + gap_penalty

            best = max(diag, up, left)
            score[i][j] = best

            if best == diag:
                back[i][j] = "diag"
            elif best == up:
                back[i][j] = "up"
            else:
                back[i][j] = "left"

    # TRACEBACK
    aligned1 = []
    aligned2 = []

    i, j = n, m

    while i > 0 or j > 0:
        direction = back[i][j]

        if direction == "diag": # diag represents a match
            aligned1.append(seq1[i - 1])
            aligned2.append(seq2[j - 1])
            i -= 1
            j -= 1

        elif direction == "up": # up represents a deletion in seq2
            aligned1.append(seq1[i - 1])
            aligned2.append('-')
            i -= 1

        elif direction == "left": # left represents an insertion in seq2
            aligned1.append('-')
            aligned2.append(seq2[j - 1])
            j -= 1

    aligned1.reverse()
    aligned2.reverse()

    final_score = score[n][m]

    return "".join(aligned1), "".join(aligned2), final_score


def local_alignment(seq1, seq2, scoring_function):
    """Local sequence alignment using the Smith-Waterman algorithm.

    Indels should be denoted with the "-" character.

    Parameters
    ----------
    seq1: str
        First sequence to be aligned.
    seq2: str
        Second sequence to be aligned.
    scoring_function: Callable

    Returns
    -------
    str
        First aligned sequence.
    str
        Second aligned sequence.
    float
        Final score of the alignment.

    Examples
    --------
    >>> local_alignment("pending itch", "unending glitch", lambda x, y: [-1, 1][x == y])
    ('ending --itch', 'ending glitch', 9.0)

    Other alignments are not possible.

    """
    raise NotImplementedError()


## This is an example scoring function, you should implement a version which uses a scoring matrix 
def scoring_function_simple(aa_i,aa_j):
    score = [-1, 1][aa_i == aa_j]
    return (score)
